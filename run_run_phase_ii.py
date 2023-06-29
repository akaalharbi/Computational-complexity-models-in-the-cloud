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



def times_required():
    """Given all attack parameters estimate how long does it take.

    This inlcude regenerating the message and collecting enough collisions.-
    """
    pass

def get_server_names():
    """Get the server names to retrieve their power usage later.

    This method will also deduce the number of servers!
    """
    import os
    with open(os.environ["OAR_NODEFILE"], "r") as f:
        server_names = f.readlines()

    server_names = [name.strip() for name in set(server_names)]

    return server_names

def copy_folder(): # todo args
    pass

def attack_choices():  # todo set parameters
    """Return the best five attack parameters on a given and a given nservers.
    """
    pass

if __name__ == "__main__":
    # get the number of servers.
    # find the best 5 attack parameters.
    # for each one of them create a folder_name_parameters
    # get the estimated time of the attack given a set of parameters
    # run each attack in its repsoecting folder for at most for at
    # most 3x the time estimated time.
    # collect the energy consumption.
    pass
