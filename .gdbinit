# command:  mpirun -np 2 emacs --eval '(gdb "gdb -i=mi -x .gdbinit ./phase_ii")'

# break sender.c:extract_dist_points

break receiver.c:103


break sender.c:540
# commands
#   call printf("idx in word_buf = %lu\n", nbytes_per_server*server_id + local_idx*one_elm_size)
#   call printf("ctr=%llu while copied ctr=%llu\n", msg_ctrs[i], ((u64*) my_M)[0])
#   call puts("let's compare digests: 1st orginal, 2nd copied")
#   call print_char(&digests[i*N], N)
#   call print_char(&work_buf[nbytes_per_server*server_id + local_idx * one_elm_size + sizeof(CTR_TYPE)], N)
  
# end



