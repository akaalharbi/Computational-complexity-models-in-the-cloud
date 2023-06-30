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
    server_names = get_server_names()
    # combine server names
    url_to_download = f"{server_names}"

    # download url
    # save the data into a file


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
    """Return the best five attack parameters on a given and a given nservers.
    """
    # loop over available choices
    # difficulty <= 8
    choices = []  # (nstates, nsenders, nreceivers, time_needed)

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
    os.system("mkdir -p experiments")

    # download the power consumption data
    # return to the base folder after the end of each experiment
    os.chdir("../../")
