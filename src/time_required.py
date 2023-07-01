"""Estimate the time required to get a collision.

Given number of senders, receivers, n, and l.
"""


def seconds_2_time(t):
    """Convert times given as seconds into a readable string."""
    from math import floor

    t = float(t)
    days = floor(t/(3600*24))
    t = t - days*24*3600
    hours = floor(t/3600)
    t = t - hours*3600
    minutes = floor(t/60)
    t = t - minutes*60

    return f"{days} days, {hours} hours, {minutes} mins, {floor(t)} sec"


def regen_msg_time(nsenders,
                   nreceivers,
                   hashes_sec_core,
                   dict_add_sec,
                   difficulty,
                   nhashes_stored):
    """Return number of seconds needed to regenerate the long message."""

    nsecs_sender = nhashes_stored / (nsenders*hashes_sec_core)
    nsecs_receiver = (nhashes_stored/(2**difficulty))
    nsecs_receiver = nsecs_receiver / (nreceivers*dict_add_sec)

    return max(nsecs_receiver, nsecs_sender)


def nqueries_sender(nsenders, hashes_sec_core, difficulty):
    """Return how many queries senders can generate per second."""
    return nsenders*hashes_sec_core/(2**difficulty)


def nqueries_receiver(nreceivers, dict_queries_sec):
    """
    Return how many queries receivers can make in a second
    """
    return nreceivers * dict_queries_sec


def nqueries(n,  # as number of bits
             nsenders,
             nreceivers,
             nhashes_stored,
             hashes_sec_core,
             dict_queries_sec,
             difficulty):
    """How many seconds needed to get reach the expected number of queries."""
    nqueries_sec = min(nqueries_sender(nsenders, hashes_sec_core, difficulty),
                       nqueries_receiver(nreceivers, dict_queries_sec))

    return 2**n / (nqueries_sec*nhashes_stored)


def times_required(n,  # n in bits not bytes
                   nstates,  # nstates stored in phase_i
                   nsenders,
                   nreceivers,
                   difficulty):
    """Given all attack parameters estimate how long does it take.

    This inlcude regenerating the message and collecting enough collisions.-
    """
    hashes_sec_core = 2**25.015134275003597
    dict_queries_sec = 2**21.863350
    # how many hashes can oure compressed file
    # hashes_sec_phase_i = 2^25  # 2^24.72
    dict_add_sec = 2**23.41
    nhashes_stored = nstates

    # todo check the arguments are in correct order!
    time_rgen_msg = regen_msg_time(nsenders,
                                   nreceivers,
                                   hashes_sec_core,
                                   dict_add_sec,
                                   difficulty,
                                   nhashes_stored,
                                   )

    time_phase_ii = nqueries(n, nsenders, nreceivers, nhashes_stored,
                             hashes_sec_core, dict_queries_sec, difficulty)

    # total time required in seconds
    return time_rgen_msg + time_phase_ii
