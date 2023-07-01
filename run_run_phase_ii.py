"""
Given a number of servers, run phase_ii using different parameters, and prepare
files to be used in phase_iii.

This program will:
- take  number of servers as input, N_min, N_max, t time
- it will pick the best 5 attack parameters, create a folder for each.
- cut the number of states and adjust the difficulty.
- download the energy consumption to the data folder.
-----------------------------------------------------
Planning for run_run_phase_iii:
* todo!
"""
from src.time_required import time_required






def get_server_names():
    """Get the server names to retrieve their power usage later.

    This method will also deduce the number of servers!
    """
    import os
    with open(os.environ["OAR_NODEFILE"], "r") as f:
        server_names = f.readlines()

    server_names = [name.strip() for name in set(server_names)]

    return server_names

def get_power_consumption_data(t_start, t_end):
    import os
    job_id = os.environ["OAR_JOBID"]

    # &start_time=2021-06-08T15:00&end_time=2021-06-08T17:00
    # combine server names
    url_to_download = f"https://api.grid5000.fr/stable/sites/grenoble/metrics?job_id={job_id}&metrics=bmc_node_power_watt&&start_time={t_start}&end_time={t_end}"

    commad = f"curl '{url_to_download}' | jq -r '.[] | [.timestamp, .device_id, .metric_id, .value, .labels|tostring] | @csv'  > data/bmc_watt{job_id}.csv"
    # download & save the data into a file
    os.system(commad)
    # curl 'https://api.grid5000.fr/stable/sites/grenoble/metrics?job_id=2296184&metrics=bmc_node_power_watt' | jq -r '.[] | [.timestamp, .device_id, .metric_id, .value, .labels|tostring] | @csv'  > bmc_watt.cs




def init_folder(n,
                nstates,
                nsenders,
                nreceivers,):
    """Copy the long message attack template folder to modify the new foldre.
    """

    import os
    os.system(f"cp -r long-message-attack experiments/{n}_nstates{nstates}_nsenders{nsenders}_nreceivers{nreceivers}")
    os.chdir(f"experiments/{n}_nstates{nstates}_nsenders{nsenders}_nreceivers{nreceivers}")

    # truncate the states file
    os.system(f"truncate --size={nstates*32} data/states")
    


def attack_choices(n,
                   nservers,
                   ncores_per_server,
                   server_memory):  # todo set parameters
    """
    Return the best five attack parameters on a given and a given nservers.
    """

    import os
    # loop over available choices
    # difficulty <= 8
    choices = []  # (nstates, nsenders, nreceivers, time_needed)
    INTERVAL = 2**23
    # how many bytes in the file
    nbytes = os.stat("data/states").st_size
    # how many states when the file gets uncompressed
    # n_available_states = (nbytes/32) * (INTERVAL)
    total_memory_nstates = nbytes*INTERVAL
    available_memory = server_memory*nservers

    print(available_memory, total_memory_nstates)

    for nreceivers in range(nservers,  # start
                            nservers*ncores_per_server - nservers,
                            nservers):  # step size
        nsenders = (ncores_per_server*nservers) - nreceivers
        for diff in range(0, 9):
            # compute l

            # rule 1: don't pass n/2 limit
            # rule 2: we can't use more states than memory allows us!
            # rule 3: we can't fill memory with nstates less than available
            #         in phase_i!
            nstates = min(2**(n/2),
                          total_memory_nstates/(2**diff),
                          available_memory/(32))

            t = time_required(n,
                              nstates,
                              nsenders,
                              nreceivers,
                              diff)

            choices.append((t, nsenders, nreceivers, nstates, diff))

    # sort choices according to the time in ascending order
    sorted(choices, key=lambda tup: tup[0])
    return choices[:5]


if __name__ == "__main__":
    import os
    # get the number of servers.
    # find the best 5 attack parameters.
    # for each one of them create a folder_name_parameters
    # get the estimated time of the attack given a set of parameters
    # run each attack in its repsoecting folder for at most for at
    # most 3x the time estimated time.
    # collect the energy consumption.
    os.system(f"python run_phase_ii.py -N {1}")
    os.system("mkdir -p experiments")

    # download the power consumption data
    # return to the base folder after the end of each experiment
    os.chdir("../../")
