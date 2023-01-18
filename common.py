import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
import mysql.connector
import bisect

def get_info_from_fname(fname):
    months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
    cell_name = ['ypd', 'water', 'Jurkat']
    t1900 = datetime.strptime("00", "%S")

    expr_info = {}
    fname_tokens = [s for s in fname.split('.')[0].split('/')[-1].split('_')]
    chip_lst = [s for s in fname_tokens if s[:4] == 'chip']
    chip_str = str(int(chip_lst[0][4:])) if len(chip_lst) == 1 else 'NULL'
    if len(chip_lst) == 1:
        expr_info['chip_id'] = str(int(chip_lst[0][4:]))

    date_token = [s for s in fname_tokens if s.lower()[:3] in months]
    date = datetime.strptime(date_token[0], "%b%d") if len(date_token) == 1 else None
    year_offset = (t1900 + (datetime.now() - date)).year
    date_new = datetime.strptime(str(year_offset) + date.strftime("/%m/%d"), "%Y/%m/%d")
    expr_info['date'] = date_new

    cell_included = [s for s in cell_name if s in fname]
    if len(cell_included) == 1:
        expr_info['cell_name'] = cell_included[0]

    return expr_info

def count2af(count, clk_freq = 4e6, sampling_cycle = 2**22, sensitivity = 590):
    time_sampled = sampling_cycle * (1 / clk_freq) # In unit second
    freq_khz = 1e-3 * count / time_sampled
    cap_ff = (1 / sensitivity) * freq_khz
    return cap_ff * 1000

def count2hz(count, clk_freq = 4e6, sampling_cycle = 2**22, sensitivity = 590):
    time_sampled = sampling_cycle * (1 / clk_freq) # In unit second
    freq_hz = count / time_sampled
    return freq_hz

def mysql_send(cursor, query, log=True):
    assert type(query) is str
    if log:
        print("> ", query)
    cursor.execute(query)
    ret = cursor.fetchall()
    if log:
        print("[mySQL]", ret)
    return ret

def find_same_expr(conditions, cursor):
    query = 'SELECT * FROM expr WHERE %s;' % conditions
    print(query)
    ret = mysql_send(cursor, query)
    expr_samecond = [r[0] for r in ret]

    data_len = len(sensor_data_df) * 25
    expr_same = None
    for idx_expr in expr_samecond:
        query = 'SELECT * FROM capval WHERE expr_id = %d;' % idx_expr
        ret = mysql_send(cursor, query, log=False)
        # print("Expr %d (%d) vs. date_len (%d)" % (idx_expr, len(ret), data_len))
        if len(ret) == data_len:
            print("Expr %d is the same" % idx_expr)
            assert expr_same is None
            expr_same = idx_expr

    return expr_same

def decsim(d0, d1):
    bin0 = bin(int(d0)).zfill(32)
    bin1 = bin(int(d1)).zfill(32)
    return sum([1 if c0 == c1 else 0 for c0, c1 in zip(bin0, bin1)])

def get_percentile(x, base):
    return sum(x > base) / len(x > base)
def normalize(x, epsilon=1e-18, x_scale=None):
    assert x_scale is None or x_scale == x_scale
    wlen = x.shape[0] // 10
    if x_scale is None:
        x_scale_rcp = np.percentile(np.diff(x), 75) - np.percentile(np.diff(x), 25) + epsilon
        x_scale = 1 / x_scale_rcp
    # x_scale = np.percentile(x, 75) - np.percentile(x, 25) + epsilon
    # x_mean = np.convolve(x, np.ones(wlen) / wlen, 'same') if x.size > 10 else x.mean()
    x_mean_middle = np.array(medianSlidingWindow(x, wlen))
    x_mean = np.pad(x_mean_middle, ((wlen - 1) // 2, wlen // 2))
    # plt.plot(x)
    # plt.plot(x_mean, color='lightgray')
    # plt.show()
    print("x_mean_middle.shape = ", x_mean_middle.shape, ", x_mean.shape = ", x_mean.shape)
    print(x_scale, np.percentile(x, 75),  np.percentile(x, 25), x.shape)
    z = (x - x_mean) * x_scale 
    print("z.pr25 = ", np.percentile(z, 25), ", z.pr75 = ", np.percentile(z, 75))
    return z

def medianSlidingWindow(nums, k):
    if k == 0: return []
    ans = []
    window = sorted(nums[0:k])
    for i in range(k, len(nums) + 1):
      ans.append((window[k // 2] + window[(k - 1) // 2]) / 2.0)
      if i == len(nums): break
      index = bisect.bisect_left(window, nums[i - k])
      window.pop(index)
      bisect.insort_left(window, nums[i])
    return ans

def unit2scale(unit, clk_freq=4e6, sampling_cycle=2**22):
    if unit == 'aF':
        scale = count2af(1, clk_freq=clk_freq, sampling_cycle=sampling_cycle)
    elif unit == 'fF':
        scale = count2af(1, clk_freq=clk_freq, sampling_cycle=sampling_cycle) * 1e-3
    elif unit.lower() == 'hz':
        scale = count2hz(1, clk_freq=clk_freq, sampling_cycle=sampling_cycle)
    elif unit.lower() == 'khz':
        scale = count2hz(1, clk_freq=clk_freq, sampling_cycle=sampling_cycle) * 1e-3
    elif unit == 'count':
        scale = 1
    else:
        raise Error()
    return scale

def median_pool(cap_dfs, keyy="deltac"):
    assert len(cap_dfs) > 0
    t_ref = cap_dfs[0]["time"]
    v_aligned = []
    for cap_df in cap_dfs:
        mask = cap_df["zscore"].abs() < 2
        t = cap_df["time"][mask]
        v = cap_df[keyy][mask]
        v_aligned.append(np.interp(t_ref, t, v))
    v2d = np.array(v_aligned)
    assert t_ref.size == v2d.shape[1]
    v2d_median = np.median(v2d, axis=0)
    assert t_ref.size == v2d_median.size
    return t_ref, v2d_median
    import pdb; pdb.set_trace()

def median_filt(data, ma_len):
    # Type check
    assert type(data) is np.ndarray or type(data) is list, print("type(data) = ", type(data))

    # Expand (roll) with padding preprocessing
    hma = ma_len // 2
    assert 2 * hma + 1 == ma_len
    data_pad = np.pad(data, (hma, hma), 'edge')
    assert len(data_pad) == len(data) + 2 * hma
    data_unroll = np.array([np.roll(data_pad, x) for x in range(-hma, hma + 1)])
    assert data_unroll.shape == (ma_len, len(data) + ma_len - 1), print("data_unroll.shape = ", data_unroll.shape, ", should be (%d, %d)" % (ma_len, len(data) + ma_len - 1))

    # Get median from expanded array
    result = np.median(data_unroll, axis=0)[hma:-hma] if hma > 0 else np.median(data_unroll, axis=0)
    assert len(result) == len(data), print("len(result) = %d, len(data) = %d" % (len(result), len(data)))
    return result

def exp_filt(data_raw, ma_len):

    # Type check
    data = data_raw.tolist()
    assert type(data) is np.ndarray or type(data) is list, print("type(data) = ", type(data))

    # Expand (roll) with padding preprocessing
    hma = (ma_len + 1) // 2
    assert 2 * hma == ma_len + 1
    hhalf = [np.exp(i) for i in range(hma - 1)]
    h_unnormalized = hhalf + [np.exp(hma-1)] + [i for i in reversed(hhalf)]
    h = [x / sum(h_unnormalized) for x in h_unnormalized]
    assert len(h) == ma_len
    assert sum(h) > 0.98 and sum(h) < 1.02, print("sum(h) = ", sum(h))

    # Get median from expanded array
    data_padded = (hma - 1) * [data[0],] + data + (hma - 1) * [data[-1],]
    result = np.convolve(data_padded, h, 'valid')
    assert len(result) == len(data), print("len(data) = %d, len(result) = %d" % (len(data), len(result)))
    return result

def get_tmask(t, tmin: float, tmax: float):
    mask = np.ones_like(t)
    if tmin is not None:
        mask[t < tmin] = 0
    if tmax is not None:
        mask[t > tmax] = 0
    return mask

