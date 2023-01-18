from common import normalize, decsim, mysql_send, count2af
def query_expr(cursor, expr_id, min_similarity=28, allow_zero=False, cmin=None, cmax=None, cond = "", scale=1, db_name='capval'):
    query = 'DESCRIBE %s;' % db_name
    query_colname = [x[0] for x in mysql_send(cursor, query, log=False)]
    print("query_colname = ", query_colname)
    idx_value = query_colname.index('value')
    idx_sid = query_colname.index('sensor_id')
    idx_time = query_colname.index('time')
    cap_value = {}
    # for i in range(64):
    #     cap_value[i] = []
    query = 'SELECT * FROM %s WHERE expr_id = %d%s;' % (db_name, expr_id, cond)
    all_expr = mysql_send(cursor, query, log=False)
    assert len(all_expr) > 0, print("Query: %s" % query)
    vprev = 0
    t0 = all_expr[0][idx_time]
    print("Total data: %d points" % len(all_expr))
    for e in all_expr:
        v = e[idx_value]
        if type(v) is int and decsim(vprev, v) <= min_similarity:
            if allow_zero or (v != 0):
                vprev = v
                v_scaled = v * scale
                sid = e[idx_sid]
                assert sid is not None, print("Sensor id is not in DB")
                assert type(sid) is int
                # if sid < 64:
                #     if cmin is None or v_scaled > cmin:
                #         if cmax is None or v_scaled < cmax:
                #             t = e[idx_time] - t0
                #             cap_value[sid].append((t, v_scaled))
                if sid in cap_value.keys():
                    t = e[idx_time] - t0
                    cap_value[sid].append((t, v_scaled))
                else:
                    print("Sensor #%d doesn't exist by default" % sid)
                    cap_value[sid] = []
                    t = e[idx_time] - t0
                    cap_value[sid].append((t, v_scaled))
                    
        elif type(v) is float:
            sid = e[idx_sid]
            t = e[idx_time] - t0
            cap_value[sid].append((t, v))
    if 23 not in cap_value:
        print("23 not exist (gen-1 chip) -> duplicate from 22")
        cap_value[23] = cap_value[22].copy()
    print("Filter removed %d zero point" % len([e for e in all_expr if e[idx_value] == 0]))
    print("After filter: %d points" % sum([len(cap_value[i]) for i in range(25)]))
    return cap_value

if __name__ == "__main__":
    import mysql.connector
    from query_desc import query_desc, query_chipid, get_chip_gen
    from common import unit2scale
    import argparse

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
    args = parser.parse_args()
    
    cnx = mysql.connector.connect(user='root-icbio', database='cell_capacitance', password='Pessw0rd-', host='10.0.0.144')
    cursor = cnx.cursor()

    for eid in args.expr_id:
        min_similarity = 28
        chip_gen = get_chip_gen(cursor, eid)
        print("chip_gen = ", chip_gen)
        if chip_gen == 1:
            clk_freq, sampling_cycle = 4e6, 2 ** 22
        elif chip_gen == 2:
            clk_freq, sampling_cycle = 1, 1
        scale = unit2scale(args.unit, clk_freq=clk_freq, sampling_cycle=sampling_cycle)
        cap_value = query_expr(cursor, eid, scale=scale, cmin=args.cmin, cmax=args.cmax)

        if 1: # First round
            to_be_dumped = []
            for k in cap_value.keys():
                values = [v for t, v in cap_value[k]]
                to_be_dumped.append(values[0])
            print("values = ", to_be_dumped)
            fname = "dumped_csv/expr%d-firstround.csv" % eid
            print("Export %s" % fname)
            with open(fname, 'w') as f:
                f.write("\n".join(["%.1f" % v for v in to_be_dumped]))
    import pdb; pdb.set_trace()
