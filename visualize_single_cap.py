import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mysql.connector
from common import normalize, decsim, mysql_send, unit2scale, median_pool, median_filt, exp_filt, get_tmask
from query_expr import query_expr
from query_desc import query_desc, query_chipid
from scipy.signal import savgol_filter, find_peaks
import matplotlib
from scipy.ndimage import gaussian_filter1d
matplotlib.use('TKAgg')

def unit2ylabel(unit):
    unit_lower = unit.lower()
    if 'hz' in unit_lower:
        return "$\Delta$ f (%s)" % unit
    elif 'af' in unit_lower or 'ff' in unit_lower:
        return "$\Delta$ C (%s)" % unit
    else:
        raise ValueError

def query_expr_specific(cursor, expr_id, name):
    query = 'DESCRIBE capval;'
    query_colname = [x[0] for x in mysql_send(cursor, query, log=False)]
    idx_value = query_colname.index(name)
    idx_sid = query_colname.index('sensor_id')
    idx_time = query_colname.index('time')
    cap_value = {}
    for i in range(25):
        cap_value[i] = []
    query = 'SELECT * FROM capval WHERE expr_id = %d;' % expr_id
    all_expr = mysql_send(cursor, query, log=False)
    t0 = all_expr[0][idx_time]
    for e in all_expr:
        v = e[idx_value]
        sid = e[idx_sid]
        t = e[idx_time] - t0
        cap_value[sid].append((t, v))
    cap_value[23] = cap_value[22].copy()
    return cap_value

def get_df_avg(cap_df, args):
    print(cap_df.keys(), len(cap_df.keys()))
    keys = [kk for kk in cap_df.keys() if type(kk) is int or kk.isnumeric()]
    if len(keys) == 64:
        chip_gen = 2
        nsensor = 64
    elif len(keys) == 25:
        chip_gen = 1
        nsensor = 16
    else:
        msg = "Cannot recognize chip with %d sensosrs" % len(keys)
        raise ValueError(msg)
    t_avg = cap_df[11]["time"]
    df = pd.DataFrame({"time": t_avg})
    for k in cap_df[0].keys():
        v_total = np.zeros_like(t_avg)
        for i in range(nsensor):
            if i not in cap_df:
                print ("Warning: Sensor #%d is not included" % i)
                continue
            mask = cap_df[i]["zscore"].abs() < args.zscore_th
            if mask.sum() == 0 and len(mask) > 0:
                import pdb; pdb.set_trace()
                raise Exception("Everything is masked")
            selected = cap_df[i][mask]
            err_msg = "sensor #%d is empty" % i
            assert len(selected["time"].tolist()) > 0, print(err_msg)
            v_total += np.interp(t_avg, selected["time"], selected[k])
        df[k] = v_total / nsensor
    return df

def get_df_temp(cap_df, t, temp, args):
    df = pd.DataFrame({"time": t, "value": temp, "value-antitemp": temp})
    if args.t_offset is not None:
        idx0 = max(sum(np.array(t) < args.t_offset), 3)
        v0 = np.median(temp[idx0-3:idx0+3])
        df['deltac'] = temp - v0
        df['deltac-antitemp'] = temp - v0
    if args.evarate is not None:
        df['value-antieva'] = df['value'] + df['time'] * args.evarate
        df['value-antieva-antitemp'] = df['value'] + df['time'] * args.evarate
        if 'deltac' in df:
            df['deltac-antieva'] = df['deltac'] + df['time'] * args.evarate
            df['deltac-antieva-antitemp'] = df['deltac'] + df['time'] * args.evarate

def get_cap_df_from_cap_value(cap_value, args):
    keyy_raw = "value" if args.t_offset is None else "deltac"
    cap_df = {}
    lmax = max([len(cap_value[k]) for k in cap_value])

    # Load all the sensors
    cap_value_np = []
    for k in cap_value:
        print(k, np.std([v for t, v in cap_value[k]]))
        cap_value_np.append([v for t, v in cap_value[k]])
    print("std = ", [np.std(v) for v in cap_value_np])
    cap_std = np.mean([np.std(v) for v in cap_value_np])
    print("cap_std = ", cap_std)
    # import pdb; pdb.set_trace()

    for k in cap_value:
        if len(cap_value[k]) < min(10, lmax // 2):
            continue
        t_all = np.array([t.total_seconds() / 3600 for t, v in cap_value[k]])
        v_all = np.array([v for t, v in cap_value[k]])
        t_min, t_max = args.tmin or min(t_all), args.tmax or max(t_all)
        tmask = (t_all >= t_min) * (t_all <= t_max)
        assert t_all[tmask].size > 0
        df = pd.DataFrame({"time": t_all[tmask], "value": v_all[tmask]})
        df['zscore'] = normalize(v_all[tmask], x_scale = 1 / cap_std)

        if args.t_offset is not None:
            idx0 = max(sum(np.array(t_all) < args.t_offset), 9)
            v0 = np.median(v_all[idx0-9:idx0+9])
            df['deltac'] = df['value'] - v0

        if args.evarate is not None:
            df['value-antieva'] = df['value'] + df['time'] * args.evarate
            if 'deltac' in df:
                df['deltac-antieva'] = df['deltac'] + df['time'] * args.evarate
        cap_df[k] = df

    # Anti-reference
    if 22 in cap_df:
        cap_df_ref = cap_df[22]
        tref = cap_df_ref['time']
        for k in cap_value:
            if k < 16: break
            tk = cap_df[k]['time']
            cap_df[k]['value-antiref'] = cap_df[k]['value'] - np.interp(tk, tref, cap_df_ref['value'])
            if args.t_offset is not None:
                cap_df[k]['deltac-antiref'] = cap_df[k]['deltac'] - np.interp(tk, tref, cap_df_ref['deltac'])
            if args.evarate is not None:
                cap_df[k]['value-antieva-antiref'] = cap_df[k]['value-antieva'] - np.interp(tk, tref, cap_df_ref['value-antieva'])
                if args.t_offset is not None:
                    cap_df[k]['deltac-antieva-antiref'] = cap_df[k]['deltac-antieva'] - np.interp(tk, tref, cap_df_ref['deltac-antieva'])

    # Add "average" entry
    if args.show_average:
        cap_df["average"] = get_df_avg(cap_df, args)

    # Add "temperature"
    if args.temperature is not None:
        assert type(args.temperature) is list
        t, temp = median_pool([cap_df[i] for i in args.temperature], keyy=keyy_raw)

        # Add temperature-related entry for other sensors
        for k in cap_value:
            tk = cap_df[k]["time"]
            vk = cap_df[k]["value"]
            tempk = np.interp(tk, t, temp)
            if args.temprate is not None:
                if args.evarate is not None:
                    cap_df[k]["value-antieva-antitemp"] = cap_df[k]["value-antieva"] - args.temprate * tempk
                cap_df[k]["value-antitemp"] = vk - args.temprate * tempk
                if args.t_offset is not None:
                    cap_df[k]["deltac-antitemp"] = cap_df[k]['deltac'] - args.temprate * tempk
                    if args.evarate is not None:
                        cap_df[k]["deltac-antieva-antitemp"] = cap_df[k]["deltac-antieva"] - args.temprate * tempk

        # Add "temperature" entry
        cap_df["temperature"] = get_df_temp(cap_df, t, temp, args)
        
    # Add "average"
    if args.show_average:
        cap_df["average"] = get_df_avg(cap_df, args)

    return cap_df

def get_chip_gen(cursor, eid):
    chip_id = query_chipid(cursor, eid)
    print("chip_id = ", chip_id)
    return 2 if chip_id > 100 else 1

sensor_fit = [[1.00, 0], [0.97, 178], [0.98, 13], [0.99, 218], [1.00, 80], [1.01, 21], [0.97, 51], [1.00, -61], [0.99, -61], [0.97, 79], [0.98, -9], [0.99, 27], [0.98, -42], [0.99, -37], [1.00, -50], [0.99, 19]]

def get_mask_without_spike(values):

    diff_raw = np.diff(values, axis=1).mean(0)
    assert diff_raw.shape == (values.shape[1] - 1,)
#     diff_raw = np.zeros_like(values[0,1:])
# #     diff_raw = np.abs(df[0][keyy][1:].to_numpy() - cap_df[0][keyy][:-1].to_numpy())
#     for i in range(values.shape[0]):
# #         y_interp = np.interp(cap_time, cap_df[i]["time"], cap_df[i][keyy])
#         diff_raw += np.diff(values)
    diff_mask = (diff_raw < np.percentile(diff_raw, 99))  & (diff_raw > np.percentile(diff_raw, 1))
    diff = diff_raw[diff_mask]
    diff_mu5sigma = diff.mean() + 5 * diff.std()
    diff_mu3sigma = diff.mean() + 3 * diff.std()
    diff_mu2sigma = diff.mean() + 2 * diff.std()
    print("diff.std() = ", diff.std())
    print("diff_mu2sigma = ", diff_mu2sigma)
    spikes = find_peaks(diff_raw, prominence=diff_mu5sigma, height=diff_mu3sigma, rel_height=.1)[0]
    spikes = []
    print("spikes = ", spikes)
#     import pdb; pdb.set_trace()
    print(find_peaks(diff_raw, prominence=diff_mu5sigma, height=diff_mu3sigma, rel_height=.1)[1])
    # plt.plot(diff_raw)
    mask_nospike = np.ones_like(values[0,:], dtype=bool)
    for spk in spikes:
        if spk > 30:
            print("spike at = ", cap_time[spk])
            mask_nospike[spk-5:spk+30] = False
    #     plt.scatter(spk, diff_raw[spk])
    # plt.show()
    return mask_nospike

def get_value_aligned(cap_df, t, keyy, indices=[0, 1, 2, 3, 4, 5]):
    assert len(t.shape) == 1
    values_lst = []
    for idx in indices:
        y_interp = np.interp(t, cap_df[idx]["time"], cap_df[idx][keyy])
        values_lst.append(y_interp)
    values = np.array(values_lst)
    assert values.shape == (len(indices), t.shape[0])
    return values


    diff_raw = np.abs(cap_df[0][keyy][1:].to_numpy() - cap_df[0][keyy][:-1].to_numpy())
    for i in range(1, 5):
        y_interp = np.interp(cap_time, cap_df[i]["time"], cap_df[i][keyy])
        diff_raw += np.abs(y_interp[1:] - y_interp[:-1])

def get_df_time(cap_df, t_resolution=0):

    # Get uniform t-axis value
    cap_time_lst = [cap_df[0]["time"].tolist()[0]]
    for tt in cap_df[0]["time"]:
        if tt - cap_time_lst[-1] > t_resolution:
            cap_time_lst.append(tt)
        # while tt - cap_time[-1] > args.plt_resolution:
        #     cap_time.append(cap_time[-1] + args.plt_resolution)
    cap_time = np.array(cap_time_lst)
    # if args.dump_csv:
    #     to_be_dumped['cap_time'] = cap_time
    # cap_time = cap_df[0]["time"].to_numpy()
    print("cap_time = ", ", ".join(["%.1f" % tt for tt in cap_time]))
    assert (cap_time[1:] - cap_time[:-1]).min() > -1e-3
    return cap_time

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--expr-id', type=int, nargs='+', default=[])
    parser.add_argument('--output', type=str, default=None)
    parser.add_argument('--title', type=str, default=None)
    parser.add_argument('--t-offset', type=int, default=None)
    parser.add_argument('--unit', type=str, default='aF')
    parser.add_argument('--average', action='store_true')
    parser.add_argument('--temperature', type=int, nargs="*", default=None)
    parser.add_argument('--allow-zero', action='store_true')
    parser.add_argument('--show-value', type=float, nargs='+', default=[])
    parser.add_argument('--grp1', type=int, nargs='+', default=[])
    parser.add_argument('--grp2', type=int, nargs='+', default=[])
    parser.add_argument('--grp3', type=int, nargs='+', default=[])
    parser.add_argument('--grp4', type=int, nargs='+', default=[])
    parser.add_argument('--grp5', type=int, nargs='+', default=[])
    parser.add_argument('--grp6', type=int, nargs='+', default=[])
    parser.add_argument('--grp7', type=int, nargs='+', default=[])
    parser.add_argument('--grp8', type=int, nargs='+', default=[])
    parser.add_argument('--grp9', type=int, nargs='+', default=[])
    parser.add_argument('--grp10', type=int, nargs='+', default=[])
    parser.add_argument('--grp11', type=int, nargs='+', default=[])
    parser.add_argument('--grp12', type=int, nargs='+', default=[])
    parser.add_argument('--grp13', type=int, nargs='+', default=[])
    parser.add_argument('--grp14', type=int, nargs='+', default=[])
    parser.add_argument('--grp15', type=int, nargs='+', default=[])
    parser.add_argument('--grp16', type=int, nargs='+', default=[])
    parser.add_argument('--show-is-water', action='store_true')
    parser.add_argument('--show-data-point', action='store_true')
    parser.add_argument('--show-average', action='store_true')
    parser.add_argument('--show-average-raw', action='store_true')
    parser.add_argument('--show-temperature', action='store_true')
    parser.add_argument('--anti-ref', action='store_true')
    parser.add_argument('--bittol', type=int, default=4)
    parser.add_argument('--evarate', type=float, default=None)
    parser.add_argument('--temprate', type=float, default=None)
    parser.add_argument('--sid', type=int, nargs='*', default=list(range(16)))
    parser.add_argument('--tmin', type=float, default=None)
    parser.add_argument('--tmax', type=float, default=None)
    parser.add_argument('--cmin', type=int, default=None)
    parser.add_argument('--cmax', type=int, default=None)
    parser.add_argument('--ylim-margin', type=float, default=0)
    parser.add_argument('--zscore-th', type=float, default=2)
    parser.add_argument('--ygrid-power', type=int, default=None)
    parser.add_argument('--plt-median', type=int, default=1)
    parser.add_argument('--plt-expfilt', type=int, default=1)
    parser.add_argument('--plt-resolution', type=float, default=20/3600)
    parser.add_argument('--plt-no-legend', action='store_true')
    parser.add_argument('--data-type', type=str, nargs='+', default=["cap"])
    parser.add_argument('--dump-csv', action='store_true')
    parser.add_argument('--colored', type=str, default=None)
    args = parser.parse_args()

    if args.temprate is not None:
        assert args.temperature is not None, print("No args.temperature set")
    min_similarity = 32 - args.bittol
    cmin, cmax = args.cmin, args.cmax
    keyy_raw = "value" if args.t_offset is None else "deltac"
    keyy = keyy_raw
    if args.evarate is not None:
        keyy += "-antieva"
    if args.temprate is not None:
        keyy += "-antitemp"
    if args.anti_ref:
        keyy += "-antiref"
    print("keyy_raw = ", keyy_raw)
    print("keyy = ", keyy)

    fig, ax = plt.subplots(1, 1)
    # cnx = mysql.connector.connect(user='chingyil', database='cell_capacitance', password='password', host='127.0.0.1')
    cnx = mysql.connector.connect(user='root-icbio', database='cell_capacitance', password='Pessw0rd-', host='10.0.0.144')
    cursor = cnx.cursor()

    pr2s, pr98s = [], []
    for eid in args.expr_id:

        # Get scale from chip generation and user-specified unit
        chip_gen = get_chip_gen(cursor, eid)
        print("chip_gen = ", chip_gen)
        if chip_gen == 1:
            clk_freq, smp_cycle = 4e6, 2 ** 22
        elif chip_gen == 2:
            clk_freq, smp_cycle = 1, 1
        scale = unit2scale(args.unit, clk_freq=clk_freq, sampling_cycle=smp_cycle)
        print("scale = ", scale)

#         if args.dump_csv:
#             to_be_dumped = {}
# 
        # Get valid data
        cap_value = query_expr(cursor, eid, scale=scale, cmin=cmin, cmax=cmax,\
                min_similarity=min_similarity)
        cap_df = get_cap_df_from_cap_value(cap_value, args)

        # Get uniform t-axis value
        cap_time = get_df_time(cap_df, t_resolution=args.plt_resolution)
        # cap_time_lst = [cap_df[0]["time"].tolist()[0]]
        # for tt in cap_df[0]["time"]:
        #     if tt - cap_time_lst[-1] > args.plt_resolution:
        #         cap_time_lst.append(tt)
        #     # while tt - cap_time[-1] > args.plt_resolution:
        #     #     cap_time.append(cap_time[-1] + args.plt_resolution)
        # cap_time = np.array(cap_time_lst)
        # # if args.dump_csv:
        # #     to_be_dumped['cap_time'] = cap_time
        # # cap_time = cap_df[0]["time"].to_numpy()
        # print("cap_time = ", ", ".join(["%.1f" % tt for tt in cap_time]))
        assert (cap_time[1:] - cap_time[:-1]).min() > -1e-3
        assert (cap_time[1:] > cap_time[:-1]).all()
        if args.dump_csv:
            df2dump = pd.DataFrame({'cap_time':cap_time})

        # Plot sensor signal
        cap_values_0to5 = get_value_aligned(cap_df, cap_time, keyy=keyy)
        mask_nospike = get_mask_without_spike(cap_values_0to5)
        # xxx
        # diff_raw = np.abs(cap_df[0][keyy][1:].to_numpy() - cap_df[0][keyy][:-1].to_numpy())
        # for i in range(1, 5):
        #     y_interp = np.interp(cap_time, cap_df[i]["time"], cap_df[i][keyy])
        #     diff_raw += np.abs(y_interp[1:] - y_interp[:-1])
        # diff_mask = (diff_raw < np.percentile(diff_raw, 99))  & (diff_raw > np.percentile(diff_raw, 1))
        # diff = diff_raw[diff_mask]
        # diff_mu5sigma = diff.mean() + 5 * diff.std()
        # diff_mu3sigma = diff.mean() + 3 * diff.std()
        # diff_mu2sigma = diff.mean() + 2 * diff.std()
        # print("diff.std() = ", diff.std())
        # print("diff_mu2sigma = ", diff_mu2sigma)
        # spikes = find_peaks(diff_raw, prominence=diff_mu5sigma, height=diff_mu3sigma, rel_height=.1)[0]
        # spikes = []
        # print("spikes = ", spikes)
#         import pdb; pdb.set_trace()
        # print(find_peaks(diff_raw, prominence=diff_mu5sigma, height=diff_mu3sigma, rel_height=.1)[1])
        # plt.plot(diff_raw)
        # mask_nospike = np.array(cap_time) > 0
        # for spk in spikes:
        #     if spk > 30:
        #         print("spike at = ", cap_time[spk])
        #         mask_nospike[spk-5:spk+30] = False
        #     plt.scatter(spk, diff_raw[spk])
        # plt.show()
        for i in args.sid:
            if i not in cap_df:
                print ("Warning: Sensor #%d is not included" % i)
                continue
            mask = cap_df[i]["zscore"].abs() < args.zscore_th
            assert mask.sum() > 0
            selected = cap_df[i][mask]
            # plt.subplot(2,1,1).plot(cap_df[i]["time"], cap_df[i][keyy].to_numpy())
            # # plt.subplot(2,1,2).plot(cap_df[i]["time"], cap_df[i]["zscore"].to_numpy())
            # # plt.subplot(2,1,1).scatter(cap_df[i]["time"], cap_df[i][keyy].to_numpy())
            # plt.subplot(2,1,2).scatter(selected["time"], selected["zscore"].to_numpy())
            # plt.title("Sensor #%d" % i)
            # plt.show()
            y = np.interp(cap_time, selected["time"], selected[keyy].to_numpy())
            y_interp = median_filt(gaussian_filter1d(y, args.plt_expfilt), args.plt_median)
            label = "#%2d" % i
            # ax.plot(cap_time[mask_nospike], y_interp[mask_nospike], label=label)
            tv_savgol = [(t, v) for t, v in zip(cap_time[mask_nospike], gaussian_filter1d(y_interp[mask_nospike], 9))]
            t_savgol = [t for t, v in tv_savgol]
            v_savgol = [v for t, v in tv_savgol]
            ax.plot(t_savgol, v_savgol, label=label)
            if args.dump_csv:
                # to_be_dumped[label] = np.interp(cap_time, cap_time[mask_nospike], y_interp[mask_nospike])
                # to_be_dumped[label] = np.interp(cap_time, t_savgol, v_savgol)
                df2dump.insert(1, label, np.interp(cap_time, t_savgol, v_savgol).tolist())

            # ax.plot(cap_time, y_interp * sensor_fit[i][0] + sensor_fit[i][1], label=label)
            pr2s.append(np.percentile(selected[keyy], 2))
            pr98s.append(np.percentile(selected[keyy], 98))
            if args.show_data_point:
                # plt.scatter(selected["time"], median_filt(exp_filt(selected[keyy].to_numpy(), args.plt_expfilt), args.plt_median), c='grey', s=8)
                ax.scatter(selected["time"], selected[keyy], c='grey', s=8)

        # Plot average signal (if nec.)
        if args.show_average:
            label = "(%d) " % eid + (query_desc(cursor, expr_id=eid) or "average")

            mask = cap_df["average"]["zscore"].abs() < args.zscore_th
            selected = cap_df["average"]
            diff_raw = np.abs(selected[keyy].to_numpy()[1:] - selected[keyy].to_numpy()[:-1])

            # Plot raw signal from exp/median filter
            y = np.interp(cap_time, selected["time"][mask], selected[keyy].to_numpy()[mask])
            y_interp = median_filt(gaussian_filter1d(y, args.plt_expfilt), args.plt_median)
            # ax.plot(cap_time, y_interp, color='lightgrey')
            if args.show_average_raw:
                ax.plot(selected["time"], selected[keyy], color='lightgray', zorder=1)
                pr2s.append(np.percentile(selected[keyy], 2))
                pr98s.append(np.percentile(selected[keyy], 98))
            for t in args.show_value:
                print(t, np.interp([t], cap_time, y_interp)[0])
                ax.text(t, np.interp([t], cap_time, y_interp)[0], "%.1f" % np.interp([t], cap_time, y_interp)[0])

            tv_savgol = [(t, v) for t, v in zip(cap_time[mask_nospike], gaussian_filter1d(y[mask_nospike], 9))]
            t_savgol = [t for t, v in tv_savgol]
            v_savgol = [v for t, v in tv_savgol]
            if args.colored is not None:
                color_avg = args.colored
            else:
                color_avg = 'C%d' % args.expr_id.index(eid) if len(args.expr_id) > 1 else "lightgray"
            ax.plot(t_savgol, v_savgol, label=label, c=color_avg, linewidth=3, zorder=3)
            if args.dump_csv:
                # to_be_dumped[label] = np.interp(cap_time, t_savgol, v_savgol)
                df2dump.insert(1, label, np.interp(cap_time, t_savgol, v_savgol).tolist())
            pr2s.append(np.percentile(v_savgol, 2))
            pr98s.append(np.percentile(v_savgol, 98))

        for idxg, grp in enumerate([args.grp1, args.grp2, args.grp3, args.grp4, args.grp5, args.grp6, args.grp7, args.grp8, args.grp9, args.grp10, args.grp11, args.grp12, args.grp13, args.grp14, args.grp15, args.grp16]):
            if len(grp) > 0:
                vgroup = np.zeros_like(cap_time)
                for i in grp:
                    mask = cap_df[i]["zscore"].abs() < args.zscore_th
                    assert mask.sum() > 0
                    selected = cap_df[i][mask]
                    y = np.interp(cap_time, selected["time"], selected[keyy].to_numpy())
                    vgroup += y
                label = "(%d) Group %d: Sensor %d - %d" % (eid, idxg + 1, min(grp), max(grp))
                tv_savgol = [(t, v) for t, v in zip(cap_time[mask_nospike], gaussian_filter1d((vgroup / len(grp))[mask_nospike], 9))]
                t_savgol = [t for t, v in tv_savgol]
                v_savgol = [v for t, v in tv_savgol]
                ax.plot(t_savgol, v_savgol, label=label, zorder=2)
                pr2s.append(np.percentile(v_savgol, 2))
                pr98s.append(np.percentile(v_savgol, 98))


        # Plot temperature signal (if nec.)
        if args.show_temperature:
            label = "temperature"
            selected = cap_df["temperature"]
            pr2s.append(np.percentile(selected[keyy], 2))
            pr98s.append(np.percentile(selected[keyy], 98))
            y = np.interp(cap_time, selected["time"], selected[keyy].to_numpy())
            y_interp = median_filt(gaussian_filter1d(y, args.plt_expfilt), args.plt_median)
            ax.plot(cap_time, y_interp, label=label)

        if args.dump_csv:
            fname = "dumped_csv/expr%d_df.csv" % eid
            df2dump.to_csv(fname, index=False)
            # fname = "dumped_csv/expr%d.csv" % eid
            # assert 'cap_time' in to_be_dumped
            # l = len(to_be_dumped['cap_time'])
            # keys = to_be_dumped.keys()
            #with open(fname, 'w') as f:
            #    f.write(",".join(keys) + '\n')
            #    for i in range(l):
            #        line = ",".join(["%.3f" % to_be_dumped[kk][i] for kk in keys])
            #        f.write(line + '\n')

                    


    title = "Expr ID = %s" % ",".join(["%d" % i for i in args.expr_id]) if args.title is None else args.title
    print("title = ", title)
    ax.set_title(title)
    if len(pr2s) > 0:
        diff = max(pr98s) - min(pr2s)
        print("pr2s = ", pr2s)
        print("pr2s = ", pr98s)
        ymin = min(pr2s) - args.ylim_margin * diff
        ymax = max(pr98s) + args.ylim_margin * diff
        ax.set_ylim(ymin, ymax) # min(pr2s) - 0.5 * diff, max(pr98s) + 0.5 * diff)
    ax.set_xlabel("Time(hr)")
    ylabel = unit2ylabel(args.unit)
    ax.set_ylabel(ylabel) #"-$\Delta$ C (%s)" % args.unit)
    if args.plt_no_legend:
        pass
    else:
        ax.legend()
    plt.tight_layout()
    plt.show()
    import pdb; pdb.set_trace()

if __name__ == '__main__':
    main()

