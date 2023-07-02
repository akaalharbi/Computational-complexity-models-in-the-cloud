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
# python3 run_run_phase_ii.py -n 72 --nservers 4 --ncores 32 -r 164000000000

import sys

from long_message_attack.src.time_required import time_required, seconds_2_time
INTERVAL = 2**23







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
                nreceivers,
                difficulty):
    """Copy the long message attack template folder to modify the new foldre.
    """

    import os

    os.system(f"cp -r long_message_attack experiments/N{n}_nstates{int(nstates)}_nsenders{nsenders}_nreceivers{nreceivers}_diff{difficulty}")
    os.chdir(f"experiments/N{n}_nstates{int(nstates)}_nsenders{nsenders}_nreceivers{nreceivers}_diff{difficulty}")
    
    # truncate the states file
    print(f"Going to truncate the states file to {nstates*32}")
    # since our  states file doesn't contain all states but rather
    # a compressed file
    nbytes_in_states_file = (nstates*32)//INTERVAL
    os.system(f"truncate --size={nbytes_in_states_file} data/states")



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

    # how many bytes in the file
    nbytes = os.stat("long_message_attack/data/states").st_size
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
            nstates = int(nstates)
            
            t = time_required(n,
                              nstates,
                              nsenders,
                              nreceivers,
                              diff)

            choices.append((t, nstates, nsenders, nreceivers, diff))

    # sort choices according to the time in ascending order
    # sorted(choices, 
    choices.sort(key=lambda tup: tup[0])
    return choices[:10]


if __name__ == "__main__":
    import os
    import argparse
    from math import log2
    from datetime import datetime

    # This the folder where all experiments will be done
    os.system("mkdir -p experiments")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-n",
                        type=int,
                        help="Number of BITS to be attacked, multiple of 8")

    parser.add_argument("--nservers",
                        type=int,
                        help="How many servers?")

    # assume all servers are the same.
    parser.add_argument("--ncores",
                        type=int,
                        help="How many cpu cores does a server have?")

    parser.add_argument("-r",
                        "--ram",
                        type=int,
                        help="what is the available memory for\
                        the dictionary per server?")
    args = parser.parse_args()

    # find the best 5 attack parameters.
    best_parameters = attack_choices(args.n,
                                     args.nservers,
                                     args.ncores,
                                     args.ram)

    print("summary of the parameters")
    print("-------------------------")
    total_expected_time = sum(tup[0] for tup in best_parameters)
    print(f"total exp time = {seconds_2_time(total_expected_time)}")
    for p in best_parameters:
        print(f"nsenders={p[2]}, nreceivers={p[3]},\
        l={log2(p[1])} difficulty={p[4]}")
        print(f"time={seconds_2_time(3*p[0])}")
        print("----------------------------------------")
    
    

    # for each one of them create a folder_name_parameters
    # run each attack in its repsoecting folder for at most for at
    # most 3x the time estimated time.
    for atck in best_parameters:
        print("=======================================")
        print("Now attacking: ")
        print(f"nsenders={p[2]}, nreceivers={p[3]},\
        nstates={log2(p[1])} difficulty={p[4]}")
        print(f"time={seconds_2_time(3*p[0])}")

        # N will be in bytes 
        init_folder(args.n//8,
                    atck[1],  # nstates
                    atck[2],  # nsenders
                    atck[3],  # nreceivers
                    atck[4],  # difficulty
                    )

        t_start = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        # run the attack, nstates is used from states file
        # nsenders is computed on the fly.
        N = args.n//8
        print(t_start)
        command = f"timeout {int(3*p[0])}s python run_phase_ii.py\
 --nservers {args.nservers}  --receivers {atck[3]} -N {N} --ram {args.ram} \
 --interval 23 -d {atck[4]}"
        print(command)
        os.system(command)
        t_end = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        print(t_end)
        print("************************************************")
        # collect the energy consumption.
        get_power_consumption_data(t_start, t_end)
        # return to the base folder
        os.chdir("../../")    
        # repeat
    print("we are done! many thanks to grid5000.fr")

    # todo list 
    # write another script to treat the collected data
    # write another script to run phase_iii based on the existing folder
