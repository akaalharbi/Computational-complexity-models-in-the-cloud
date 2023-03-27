# command:  mpirun -np 2 emacs --eval '(gdb "gdb -i=mi -x .gdbinit ./phase_ii")'

# break sender.c:extract_dist_points

break receiver.c:103
break sender.c:581
commands
call print_char(snd_buf, 18)
n
end
