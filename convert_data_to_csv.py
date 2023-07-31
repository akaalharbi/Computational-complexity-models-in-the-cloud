"""Convert data written by sender and receivers to csv files.

Plan:
1) Put all senders   data files into a single csv file.
2) Put all receivers data files into a single csv file.
3) Create a summary csv file that broken by phases which includes the total
   energy, total time, number of candidate, average throughput, etc.
4) Combine all summaries csv into a single csv file.
5) Check that if a folder has been processed then it should not be processed.
______________________________________________________________
# Scattered ideas:

Load balancing? how good it's?
"""

# add a function recv_msg_stat_parse() to parse a single message's
# statistics. It will return a dictionary that has all the properties.
# I am thinking of embedding this inside parse_receiver
# add a function recv_cnd_stat that parse the information about a found
# candidate.
# These two functions should write directly to the csv file?
# New thought, we should not deep nest the csv file, only parse sender writes
# to it. So, the above two functions should return a one line string to be
# added to the common csv file.


def parse_sender():
    """Process all sender_* files and store the result in sender.csv.

    This function doesn't change path, it assumes that we are in stats/ folder.
    """
    # with open(file_name, "r"): as f:
    #
    #     # else if it's a candidate 
    #     pass


def parse_receiver():
    """Process all receiver_* files and store the result in sender.csv.

    This function doesn't change path, it assumes that we are in stats/ folder.
    """
    bracket = "<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-"
    file_names = [] # get all files names that start with
    for f in file_names:
        f = open(f, "r")
        bracket_open = 0 # false
        data_extracted = ""
        for line in f:
            if line == bracket:
                # if the bracket is close then open, otherwise close it!.
                bracket_open = (bracket_open + 1)%2
            if bracket_open:
                pass  # call recv_msg_stat_parse()
            elif not bracket_open:
                pass  # call recv_cnd_stat()
            # write extracted data to the csv file


    pass


def total_consumption():
    """Get the total energy consumption and the average energy per bit.

    This function doesn't change path, it assumes that we are in stats/ folder.
    """
    pass


def create_summary():
    """Combine the 3 csv files sender, receiver, and total_consumption."""
    pass
