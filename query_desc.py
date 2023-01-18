from common import mysql_send
def query_desc(cursor, expr_id=None):
    query = 'SELECT description from expr WHERE id = %d' % expr_id
    responses = mysql_send(cursor, query, log=False)
    assert len(responses) == 1, print("len(responses) = ", len(responses))
    assert len(responses[0]) == 1, print("len(resp[0]) = ", len(responses[0]))
    return responses[0][0]

def query_chipid(cursor, expr_id):
    query = 'SELECT chip_id from expr WHERE id = %d' % expr_id
    responses = mysql_send(cursor, query, log=False)
    assert len(responses) == 1, print("len(responses) = ", len(responses))
    assert len(responses[0]) == 1, print("len(resp[0]) = ", len(responses[0]))
    return responses[0][0]

def get_chip_gen(cursor, eid):
    chip_id = query_chipid(cursor, eid)
    print("chip_id = ", chip_id)
    return 2 if chip_id > 100 else 1


