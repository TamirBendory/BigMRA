function stats = statsfun_bispect(problem, x, stats)

if mod(stats.iter, 20) == 0
    save('curr_soln_alms_TR.mat','x')
end