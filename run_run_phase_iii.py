"""Explore the runned experiments given a folder as an input."""


# 0- init folder by: creating archive, copying src & Makefile, and compiling lib folder (do we need this?)
# 1- For each folder inside the experiment folder
#   check if phase_iii has been completed, if yes skip.
#   If not, estimate how long does it take to estimate the time
#       (get number of cores, assume each core do 2^25 hash/sec)
# print the sum of the estimated time.
# 2- run phase_iii for each folder that is not completed
# 3- download energy consumption

def init_folders():
    """Copy the extra c files to the existing folder"""
    pass


