stroop <- data.table::fread(here::here("in", "ub55_stroop_behav-events.csv"))
axcpt <- data.table::fread(here::here("in", "ub55_axcpt__behav-events.csv"))
cuedts <- data.table::fread(here::here("in", "ub55_cuedts_behav-events.csv"))
stern <- data.table::fread(here::here("in", "ub55_stern_behav-events.csv"))
# subjs <- fread(here("in", "ub55_subjects.txt"))[[1]]