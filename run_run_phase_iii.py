"""Explore the runned experiments given a folder as an input."""


# 0- init folder by: creating archive, copying src & Makefile, and compiling lib folder (do we need this?)
# 1- For each folder inside the experiment folder
#   check if phase_iii has been completed, if yes skip.
#   If not, estimate how long does it take to estimate the time
#       (get number of cores, assume each core do 2^25 hash/sec)
# print the sum of the estimated time.
# 2- run phase_iii for each folder that is not completed
# 3- download energy consumption
INTERVAL = 2**13

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



def folders_that_need_treatment():
    """Check folders that have not experienced phase_III yet."""
    import os
    unworked_folders = []
    path = "experiments"
    directories = [d for d in os.listdir(path)
                   if os.path.isdir(os.path.join(path, d))]

    for d in directories:
        data_path = os.path.join(path, d) + "/data/results"

        if not os.path.exists(os.path.join(path, d) + "/data/phase_ii_status"):
            continue  # skip this folder it has not been completed by phase ii

        if not os.path.exists(data_path):
            unworked_folders.append(os.path.join(path, d))
            continue
        
        with open(data_path, "r") as f:
            for line in f:
                last_line = line
        if last_line[:3] == "END":
            # this folder passed phase_iii
            continue
        
        unworked_folders.append(os.path.join(path, d))

        # open d/data/results
        # if last line == END

    return unworked_folders


def complete_folder(path):
    """Add missing files for phase_iii to a folder."""
    import os
    os.system(f"rsync -avzP src/ {os.path.join(path, 'src/')}")
    os.system(f"rsync -avzP Makefile {os.path.join(path, 'Makefile')}")


def estimate_time_needed(path):
    """Return estimation of nseconds needed to complete phase_iii."""
    import os
    path_to_archive = os.path.join(path, "data/states")
    nbytes = os.stat(path_to_archive).st_size
    nbytes = nbytes*INTERVAL
    print(f"|{path_to_archive}|={nbytes} bytes")
    # how many seconds are estimated to regenerate the long message which is
    # the heaviest part in phase_iii
    #print((nbytes/(32*os.cpu_count()*(2**23))))
    return (nbytes/(os.cpu_count()*(2**23)))


def create_archive(path):
    """Create archive file in data/messages/."""
    import os
    path_to_messages = os.path.join(path, "data/messages/")
    print(path)
    os.system(f"rm -f {path_to_messages}archive")
    os.system(f"cat {path_to_messages}* > {path_to_messages+'archive'}")

def save_power_consumption_data(t_start, t_end):
    """Get the power consumption for the current job and save it in a file"""
    import os
    job_id = os.environ["OAR_JOBID"]

    # &start_time=2021-06-08T15:00&end_time=2021-06-08T17:00
    # combine server names
    url_to_download = f"https://api.grid5000.fr/stable/sites/grenoble/metrics?job_id={job_id}&metrics=bmc_node_power_watt&&start_time={t_start}&end_time={t_end}"

    commad = f"curl '{url_to_download}' | jq -r '.[] | [.timestamp, .device_id, .metric_id, .value, .labels|tostring] | @csv'  > data/phase_iii_bmc_watt{job_id}.csv"
    # download & save the data into a file
    os.system(commad)
    # curl 'https://api.grid5000.fr/stable/sites/grenoble/metrics?job_id=2296184&metrics=bmc_node_power_watt' | jq -r '.[] | [.timestamp, .device_id, .metric_id, .value, .labels|tostring] | @csv'  > bmc_watt.cs


if __name__ == "__main__":
    import os
    from datetime import datetime

    # relative paths: experiments/folder_name
    folders = folders_that_need_treatment()

    for d in folders:
        create_archive(d)

    time_needed = [estimate_time_needed(d) for d in folders]
    current_path = os.getcwd()

    # Printing summary information at the beginning
    print(f"There are {len(folders)} we're going to work on")
    print(f"Treating all of them is expected to take {seconds_2_time(sum(time_needed))}")

    for i in range(len(folders)):
        print(f"folder={folders[i]}, takes={time_needed[i]}sec")
    
    # Actual computation:
    for d in folders:
        complete_folder(d)


        os.chdir(d)

        # compliation
        os.system("sudo-g5k apt install nasm")
        os.system("cd lib/sha256_intel_avx/\
        && make clean\
        && make all && cd ../../")

        os.system("make clean && make all")

        # actual computation:
        t_start = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        os.system("./phase_iii")
        t_end = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

        save_power_consumption_data(t_start, t_end)
        os.chdir(current_path)  # go back to the old path

