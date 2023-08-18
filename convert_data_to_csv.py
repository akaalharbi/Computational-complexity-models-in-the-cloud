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
# WARNING: @todo senders and receivers don't agree on INTERVAL!
# see the folder long_message_attack/experiments/N10_nstates1099511627776_nsenders319_nreceivers609_diff3/data/stats$ 


def extract_receiving_stats(text):
    """Extract the data from stats of message received into a csv line."""
    import re
    # Extract the data from the text using regular expressions
    matches = re.findall(r'total=(\d+\.\d+)sec', text)
    total = matches[0]

    matches = re.findall(r'mpi_recv=(\d+\.\d+)sec', text)
    mpi_recv = matches[0]

    matches = re.findall(r'dict_add=(\d+\.\d+)sec', text)
    dict_add = matches[0]

    matches = re.findall(r'dict_add=(\d+\.\d+)sec≈2\^(\d+\.\d+)≈(\d+\.\d+)MB/sec', text)
    dict_add_speed = matches[0][2]

    matches = re.findall(r'mpi_recv=(\d+\.\d+)%', text)
    mpi_recv_percent = matches[0]

    matches = re.findall(r'dict_add=(\d+\.\d+)%', text)
    dict_add_percent = matches[0]

    matches = re.findall(r'RECV (\d+\.\d+)MB/sec', text)
    recv_speed = matches[0]

    matches = re.findall(r'exp\[all receivers\] = (\d+\.\d+) MB/sec', text)
    exp_all_receivers = matches[0]

    matches = re.findall(r'nsenders=(\d+)', text)
    nsenders = matches[0]

    matches = re.findall(r'nservers=(\d+)', text)
    nservers = matches[0]

    matches = re.findall(r'DIFFICULTY=(\d+)', text)
    difficulty = matches[0]

    matches = re.findall(r'INTERVAL=(\d+)', text)
    interval = matches[0]

    matches = re.findall(r'nmsgs_recv=(\d+)', text)
    nmsgs_recv = matches[0]

    # Combine the extracted data into a CSV line
    csv_line = f"{total},{mpi_recv},{dict_add},{dict_add_speed},{mpi_recv_percent},{dict_add_percent},{recv_speed},{exp_all_receivers},{nsenders},{nservers},{difficulty},{interval},{nmsgs_recv}"

    return csv_line


def extract_candidates_stats(text):
    """Extract number of newly found candidates and accumulated ncandidates."""
    import re

    # Extract the data from the text using regular expressions
    matches = re.findall(r'nfound_cnd=(\d+)', text)
    nfound_cnd = matches[0]

    matches = re.findall(r'new_cnd=(\d+)', text)
    new_cnd = matches[0]

    # matches = re.findall(r't=(\d+\.\d+)sec', text)
    # t = matches[0]

    # Combine the extracted data into a CSV line
    csv_line = f"{nfound_cnd},{new_cnd}"

    return csv_line


def parse_receiver_file(f_inp, f_csv):
    """Extract data from f_inp and write the extracted data to csv_name."""
    import re
    matches = re.findall(r"(\d+)", f_inp.name)
    # print(f"file_name={f_inp.name}")
    receiver_name = matches[0]
    # print(f"receiver_name={receiver_name}")


    bracket = "<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-\n"

    # read the first 6 lines which corresponds to rehashing a message.
    is_bracket_open = False
    text = ""

    # extract data for rehashing the message

    for _ in range(6):
        text += f_inp.readline()

    csv_line = extract_receiving_stats(text)
    # two entries: nnew_candidates, nfound_candidates, we ignore time
    csv_line = csv_line + ",0,0," + receiver_name + "\n"
    f_csv.write(csv_line)

    # reset variables
    text = ""
    is_bracket_open = False
    next_line_candidates = False

    print(f"DONE with regen stats {receiver_name} receiver")

    # Actual parsing:
    for line in f_inp:  # complete where stopped before
        # print(f"bracket={bracket == line}, {line}")

        if line == bracket:
            tmp = is_bracket_open
            is_bracket_open = (is_bracket_open + 1) % 2
            # print(f"is_bracket_open={is_bracket_open}")

            # we just closed a bracket that was open, i.e. receiving stats
            if tmp and not is_bracket_open:
                # treat the text
                text += line
                # print(f"text={text}")
                csv_line = extract_receiving_stats(text)
                # print("+"*40)
                # print(f"inside brackets:\n{text}")
                # print(csv_line)
                # print("+"*40)
                f_csv.write(csv_line)
                text = ""
                next_line_candidates = True

                continue  # don't do extra computation

        # or we only read a cadnidates statistics
        if line != bracket and next_line_candidates:
            text = line
            csv_line = "," + extract_candidates_stats(text)
            csv_line += f",{receiver_name}\n"
            # print("+"*40)
            # print(f"candidates {text}")
            # print(csv_line)
            # print("+"*40)
            f_csv.write(csv_line)
            # @todo add receiver name!
            text = ""  # reset the text
            next_line_candidates = False
            continue  # don't do extra computation

        # Accumlate receiving statistics text
        text += line


def parse_receivers():
    """Process all receiver_* files and store the result in sender.csv.

    This function doesn't change path, it assumes that we are in stats/ folder.
    """
    import os

    file_names = os.listdir("data/stats/")  # get all files names that start with
    file_names = filter(lambda x: "receiver_" in x, file_names)
    file_names = [os.path.join("data/stats/", f_name) for f_name in file_names]
    # print(f"files are {file_names}")
    csv_file = open("data/stats/receivers.csv", "w")
    csv_file.write("time_sec,mpi_recv_sec,dict_add_sec,dict_add_speed_MB,mpi_recv_percent,dict_add_percent,recv_speed_MB,exp_all_receivers_MB,nsenders,nreceivers,difficulty,interval,nmsgs_recv,nfound_cnd,new_cnd,receiver_name\n")
    # add csv header
    # todo

    for f in file_names:
        receiver_file = open(f, "r")
        parse_receiver_file(receiver_file, csv_file)






def extract_sender_stat(text):
    """Extract sender stats into csv format."""
    import re

    # Remove delimiters
    text = text.replace("->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->", "")

    # Split text by lines
    lines = text.strip().split('\n')

    # Initialize list to store extracted data
    data = []

    # Extract data from each line
    for line in lines:
        # Use regular expressions to extract numeric values
        values = re.findall(r'[\d.]+', line)
        data.extend(values)

    # Convert list to CSV format
    csv_data = ",".join(data)

    return csv_data


def parse_sender_file(f_inp, f_csv):
    """Loop over all senders stats and collect them inside single csv file."""
    import re
    matches = re.findall(r"(\d+)", f_inp.name)
    # print(f"file_name={f_inp.name}")
    sender_name = matches[0]
    # print(f"i.e. name={sender_name}")
    bracket = "->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->\n"
    is_bracket_open = False
    text = ""

    # The actual parsing
    for line in f_inp:
        if line == bracket:
            tmp = is_bracket_open
            is_bracket_open = (is_bracket_open + 1) % 2

            # we just closed a bracket that was open
            if tmp and not is_bracket_open:

                text += line
                # print("+++++++++++++++++++++++++++++++++")
                # print(f"{text}")
                # print("++++++++++++++++++++++++++++++++")
                csv_line = extract_sender_stat(text)
                csv_line += f",{sender_name}\n"
                # print(f"sender {sender_name} going to write\n{csv_line}")
                f_csv.write(csv_line)
                text = ""
                continue  # don't do any other computation

        text += line


def parse_senders():
    """
    Process all sender_* files and store the result in sender.csv.

    This function doesn't change path, it assumes that we are in stats/ folder.
    """
    import os

    file_names = os.listdir("data/stats/")  # get all files names that start with
    file_names = filter(lambda x: "sender_" in x, file_names)
    file_names = [os.path.join("data/stats/", f_name) for f_name in file_names]
    # print(f"files are {file_names}")
    csv_file = open("data/stats/senders.csv", "w")
    csv_file.write("time,mpi_wait_sec,hash_sec,hash_speed,hash_speed_MB,find_dist_sec,mpi_send_percent,hash_percent,find_dist_percent,send_speed_MB,exp_all_senders_speed_MB,nsenders,nreceivers,difficulty,interval,nsends,sender_name\n")
    # add csv header
    # todo

    for f in file_names:
        print(f"Going to treat {f}")
        sender_file = open(f, "r")
        parse_sender_file(sender_file, csv_file)


parse_receivers()
parse_senders()
