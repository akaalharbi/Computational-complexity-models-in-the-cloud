import os


file_path =  "statistics_parallel/37_3_29_2_stats.txt"
with open(file_path, "r") as f:

     file_log = open("log/collision_checks", "w")
     file_log.write("n l\n")
     file_log.close()

     # skip the first line
     f.readline()
     for line in f.readlines():
        line_split = line.split(",")
        n, l, idx = int(line_split[0]), int(line_split[1]), int(line_split[-2])
        path = f"messages/{n}_{l}"
        file_log = open("log/collision_checks", "a")
        file_log.write(f"n={n} l={l}, collides=")
        file_log.close()
        os.system(f"./verify_hash {path} {n} {idx}")
 
