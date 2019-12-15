function h = V_hist(velocities, centers, binwidth)

h = hist(velocities, centers);

h = h/(binwidth*sum(h));