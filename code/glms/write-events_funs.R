## about ----
## 
## functions for writing event files for afni GLMs
## 
## mike freund, 2019-10-10


onsets4afni <- function(x) {
  ## takes a vector of onsets and puts them in format for afni
  x <- x[!is.na(x)]
  if (!is.numeric(x)) stop("needs numeric!")
  if (length(x) < 1) {
    x <- "*\n"
  } else {
    x <- paste0(paste(x, collapse = " "), " \n")
  }
  return(x)
}

onsets4afni.dm <- function(x) {
  ## analogous to onsets4afni(), but adapted for use with dmBLOCK().
  if (class(x) != "data.frame") stop("needs data.frame!")
  if (any(!complete.cases(x$dm))) stop("incomplete cases in dm!")
  x <- x[complete.cases(x$time.onset), ]
  if (nrow(x) < 1) {
    ## no valid times needs -1 option, e.g.:
    ## https://afni.nimh.nih.gov/afni/community/board/read.php?1,84354,85947#msg-85947
    x <- "-1:0\n"  
  } else {
    x <- apply(x, 1, paste, collapse = ":") %>% paste0(collapse = " ") %>% paste0(" \n")
  }
  return(x)
}

onsets2file <- function(.onsets, .fname) {
  .fname <- paste0(.fname, ".txt")
  .fout <- file(.fname, "wt")
  cat(.onsets, file = .fout)
  close(.fout)
  unlink(.fout)
}



write.event <- function(
  d,
  dir.analysis,
  var.level,
  name.var,
  name.onset,
  fname.suffix,
  by.run,
  subdir.input.data = "INPUT_DATA",
  write.empty = TRUE
  ) {
  
  ## FOR DEBUGGING:
  # d <- axcpt %>% filter(subj == "132017", session == "baseline")
  # dir.analysis <- "/data/nil-external/ccp/freund/sub-subj-glms/runwise_new"
  # name.var <- "trial.type.correct"
  # var.level <- "Ang"
  # name.onset <- "time.cue.onset"
  # fname.suffix <- "correct"
  # subdir.input.data = "INPUT_DATA"
  # by.run = TRUE
  
  ## ABOUT:
  ##
  ## extracts elements of d[[colname.var]] == level to include in stimtimes file.
  ## 
  ## ARGS:
  ## 
  ## d: a dataframe. each row corresponds to a single event (trial).
  ## d only contains events for a single subject * session.
  ## must contain cols (with matching names) "subj", "session", "task", "run"
  ## 
  ## dir.analysis: the parent directory to which events are written.
  ## 
  ## var.level: a character string of length 1 that specifies the particular event type to be written.
  ## 
  ## name.var: the name of the column in d that contains elements equal to var.level, which specify the rows (events)
  ##  that should be included in the event file.
  ##
  ## name.onset: the name of the column in d that contains the onset times (actual numeric values to be written).
  ## 
  ## fname.suffix: a string to append to the end of the event file name
  ## 
  ## by.run: boolean, write separate events for each scanning run?
  ## 
  ## subdir.input.data: the subdirectory within dir.analysis that contains the subject tree of event files. default = 
  ##  "INPUT_DATA".
  ##  
  ## write.empty: boolean, if sum(d[[colname.var]] == level) == 0, should an 'empty' stimtime file be written?
  ##  Defaults to TRUE, as having an identical set of stimtime files for each subject is convenient for writing 
  ##  subsequent analysis code, as each subject will have an identically structured afni stats file.

  
  
  ## input validation
  
  d <- as.data.frame(d)
  
  names.expected <- c(name.var, name.onset, "subj", "session", "task", "run")
  names.ok <- all(names.expected %in% names(d))
  if (!names.ok) stop("d missing cols: name.var | name.onset | subj | session | task")
  
  subj <- unique(d$subj)
  sess <- unique(d$session)
  task <- unique(d$task)
  
  if (length(c(subj, sess, task, var.level)) != 4) stop(">1 subj | session | task | var.level")
  if (!sess %in% c("baseline", "proactive", "reactive")) stop("session name error")
  if (!task %in% c("Axcpt", "Cuedts", "Stern", "Stroop")) stop("task name error")
  
  ## create path / directory
  
  dir.input <- file.path(dir.analysis, subj, subdir.input.data, task, sess)
  if (!dir.exists(dir.input)) dir.create(dir.input, recursive = TRUE)
  
  ## extract relevant rows (events)
  
  is.included <- d[[name.var]] == var.level
  
  onsets1  <- d[d$run == 1 & is.included, name.onset] %>% sort %>% onsets4afni
  onsets2  <- d[d$run == 2 & is.included, name.onset] %>% sort %>% onsets4afni
  
  fname <- file.path(dir.input, paste0(subj, "_", task, "_", sess, "_", var.level, fname.suffix))
  
  if (by.run) {
    
    onsets1 %>% onsets2file(.fname = paste0(fname, "_run1"))
    onsets2 %>% onsets2file(.fname = paste0(fname, "_run2"))
    
  } else {
    
    onsets <- paste0(onsets1, onsets2, collapse = "")
    onsets %>% onsets2file(.fname = fname)
    
  }
  
  c(sum(is.included, na.rm = TRUE), fname)  ## number of onsets written, and filename (for book-keeping)

}


write.events <- function(.d, .args) {
  
  ## ABOUT
  ## a wrapper for write.event()
  ## takes a data.frame of events for a subject*task*session, and a data.frame of 
  ## values that specify arguments for write.event().
  ## writes events.
  
  subj <- unique(.d$subj)
  sess <- unique(.d$session)
  task <- unique(.d$task)
  
  if (length(c(subj, sess, task)) != 3) stop(">1 subj | session | task")
  if (!sess %in% c("baseline", "proactive", "reactive")) stop("session name error")
  if (!task %in% c("Axcpt", "Cuedts", "Stern", "Stroop")) stop("task name error")
  
  is.bad.name <- !names(.args) %in% c("var.level", "name.var", "fname.suffix", "dir.analysis", "name.onset", "by.run")
  if (any(is.bad.name)) stop(paste0("bad name(s) in .args: ", names(.args)[is.bad.name]))
  
  ## add task*session*subj info to .args
  
  .args$session <- sess
  .args$task <- task
  .args$subj <- subj
  
  ## initialize columns for loop
  
  n.events <- nrow(.args)
  .args$n.events <- as.numeric(NA)
  .args$fname <- ""
  
  for (event.i in seq_len(n.events)) {
    # event.i = 1
    
    .args[event.i, c("n.events", "fname")] <- write.event(
      d = .d,
      var.level    = .args[event.i, "var.level"],
      name.var     = .args[event.i, "name.var"],
      fname.suffix = .args[event.i, "fname.suffix"],
      dir.analysis = .args[event.i, "dir.analysis"],
      name.onset   = .args[event.i, "name.onset"],
      by.run       = .args[event.i, "by.run"]
    )
    
  }
  
  .args
  
}



write.blocks <- function(
  d,
  dir.analysis,
  by.run,
  subdir.input.data = "INPUT_DATA"
) {
  
  ## input validation
  
  d <- as.data.frame(d)
  
  names.expected <- c("subj", "session", "task", "run")
  names.ok <- all(names.expected %in% names(d))
  if (!names.ok) stop("d missing cols: subj | session | task")
  
  subj <- unique(d$subj)
  sess <- unique(d$session)
  task <- unique(d$task)
  
  if (length(c(subj, sess, task)) != 3) stop(">1 subj | session | task")
  if (!sess %in% c("baseline", "proactive", "reactive")) stop("session name error")
  if (!task %in% c("Axcpt", "Cuedts", "Stern", "Stroop")) stop("task name error")
  
  ## create directory
  
  dir.input <- file.path(dir.analysis, subj, subdir.input.data, task, sess)
  if (!dir.exists(dir.input)) dir.create(dir.input, recursive = TRUE)
  
  ## get block events
  
  block.events <- d %>% select(run, contains("time.block")) %>% filter(!duplicated(.))
  if (any(is.na(block.events))) stop("missing sustained or transient time")
  if (nrow(block.events) != 2) stop("nrow(block.events) != 2")
  
  ## bring to long:
  block.events.long <- block.events %>%
    reshape2::melt(id = "run") %>%
    mutate(
      block    = gsub(".*([1-3]).*", "\\1", variable),  ## pull out block
      variable = gsub("time.block[1-3].", "", variable)  ## pull out on / off
    ) %>% 
    tidyr::spread(variable, value) %>%
    mutate(
      dm = off - on
      ## for blocks with no off value (e.g. truncated run), replace duration with mean:
      # dm = ifelse(is.na(off) & !is.na(on), mean(dm, na.rm = TRUE), dm)
    ) %>%
    rename(time.onset = on, time.offset = off)
  
  ## build regressors
  sustained1 <- block.events.long %>% filter(run == "1") %>% arrange(time.onset) %>% select(time.onset, dm)
  sustained2 <- block.events.long %>% filter(run == "2") %>% arrange(time.onset) %>% select(time.onset, dm)
  transient1 <- block.events.long %>% filter(run == "1") %>% select(time.onset, time.offset) %>% unlist %>% sort
  transient2 <- block.events.long %>% filter(run == "2") %>% select(time.onset, time.offset) %>% unlist %>% sort
  
  fname <- file.path(dir.input, paste0(subj, "_", task, "_", sess))
  
  ## write
  
  if (by.run) {
    
    sustained1 %>% onsets4afni.dm %>% onsets2file(.fname = paste0(fname, "_block_run1"))
    sustained2 %>% onsets4afni.dm %>% onsets2file(.fname = paste0(fname, "_block_run2"))
    transient1 %>% onsets2file(.fname = paste0(fname, "_blockONandOFF_run1"))
    transient2 %>% onsets2file(.fname = paste0(fname, "_blockONandOFF_run2"))
    
  } else {
    
    sustained <- paste0(onsets4afni.dm(sustained1), onsets4afni.dm(sustained2), collapse = "")
    transient <- paste0(onsets4afni(transient1 %>% sort), onsets4afni(transient2 %>% sort), collapse = "")
    sustained %>% onsets2file(.fname = paste0(fname, "_block"))
    transient %>% onsets2file(.fname = paste0(fname, "_blockONandOFF"))
    
  }
  
}



split.movregs <- function(
  to.split,
  dir.to,
  dir.from  = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
) {
  
  ## reads movregs files from nil-bluearc, splits by run
  ## NB: THIS FUNCTION ASSUMES THE MOVREGS FILE HAS nrow == TR
  ## NB: ONLY WORKS FOR MB4 PEOPLE!
  ## 
  ## arg to.split should contain a data.frame of all subject*task*session files to split and copy.
  ## e.g., 
  # to.split <- expand.grid(
  #   subj = "132017",
  #   task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  #   session = c("baseline", "proactive", "reactive")
  #   )
  # dir.from = dl$nil.dmcc2.hcp.afni
  # dir.to   = file.path(dl$nil.external.freund, "AFNI_ANALYSIS_SUBSUBJECT")
  
  
  dir.from   <- file.path(dir.from, to.split$subj, "INPUT_DATA", to.split$task, to.split$session)
  dir.to     <- file.path(dir.to, to.split$subj, "INPUT_DATA", to.split$task, to.split$session)
  filenames  <- paste0("motion_demean_", to.split$session, ".1D")
  files.from <- file.path(dir.from, filenames)
  files.to   <- file.path(dir.to, filenames)
  
  ## create dirs if they do not exist
  for (file.i in seq_along(dir.to)) if (!dir.exists(dir.to[file.i])) dir.create(dir.to[file.i], recursive = TRUE)
  is.missing.dir <- !(dir.exists(dir.from) & dir.exists(dir.to))
  
  ## get number of TRs (for file splitting)
  trs <- expand.grid(
    task    = c("Axcpt", "Cuedts", "Stern", "Stroop"),
    session = c("proactive", "reactive", "baseline")
  )
  trs <- trs[with(trs, order(task, session)), ]
  trs$num <- c(1220, 1220, 1220, 1300, 1300, 1300, 1200, 1200, 1200, 1080, 1180, 1080)  ## order corresponds!
  
  ## read files
  has.unexpected.nrow <- vector("logical", nrow(to.split))
  for (file.i in seq_along(filenames)) {
    if (is.missing.dir[file.i]) next
    
    session.i <- gsub(".*(baseline).*|.*(proactive).*|.*(reactive).*", "\\1\\2\\3", files.from[file.i])
    task.i    <- gsub(".*(Axcpt).*|.*(Cuedts).*|.*(Stern).*|.*(Stroop).*", "\\1\\2\\3\\4", files.from[file.i])
    num.trs.i <- trs$num[trs$task == task.i & trs$session == session.i]
    
    movregs.i <- data.table::fread(
      files.from[file.i], 
      sep = " ", header = FALSE, colClasses = rep("numeric", 6), data.table = FALSE
    )
    mask.i <- data.table::fread(
      file.path(dir.from, "movregs_FD_mask.txt")[file.i],
      sep = " ", header = FALSE, colClasses = "integer", data.table = FALSE
    )[[1]]  ## all FD_mask files have same name
    
    ## remove NA values (if present) and check rownums
    movregs.i <- movregs.i[complete.cases(movregs.i), ]
    mask.i    <- mask.i[!is.na(mask.i)]
    is.unexpected.nrow <- nrow(movregs.i) != num.trs.i | length(mask.i) != num.trs.i
    if (is.unexpected.nrow) {
      has.unexpected.nrow[file.i] <- TRUE
      next
    }
    
    ## split and write
    movregs.i.run1 <- movregs.i[seq(1, num.trs.i / 2), ]
    movregs.i.run2 <- movregs.i[seq(num.trs.i / 2 + 1, num.trs.i), ]
    
    new.filename.i <- gsub("\\.1D", "", files.to[file.i])
    movregs.i.run1 %>% data.table::fwrite(paste0(new.filename.i, "_run1.1D"), sep = " ", col.names = FALSE)
    movregs.i.run2 %>% data.table::fwrite(paste0(new.filename.i, "_run2.1D"), sep = " ", col.names = FALSE)
    
    mask.i.run1 <- mask.i[seq(1, num.trs.i / 2)]
    mask.i.run2 <- mask.i[seq(num.trs.i / 2 + 1, num.trs.i)]
    
    mask.i.run1 %>% 
      data.frame %>% 
      data.table::fwrite(file.path(dir.to, "movregs_FD_mask_run1.txt")[file.i], sep = " ", col.names = FALSE)
    mask.i.run2 %>% 
      data.frame %>% 
      data.table::fwrite(file.path(dir.to, "movregs_FD_mask_run2.txt")[file.i], sep = " ", col.names = FALSE)
    
  }
  
  data.frame(to.split, is.missing.dir, has.unexpected.nrow)
  
}

