# variables for all sessions
session.ids <- c("baseline", "proactive", "reactive");
sess.ids <- c("Bas", "Pro", "Rea");  # shorter version of the session.ids; same order as session.ids
run.ids <- c(1, 2);   # AP and PA runs for each task and session

# variables for Stroop
con.ids.in <- c("Congruent", "Incongruent");   # as coded in $Congruency
con.ids.out <- c("Con", "InCon");
con.ids.cts <- rbind(c(35,19), c(35,19));  # number of trials of each type in each run, in order con.ids; 1st row run 1, 2nd row run 2

do.Stroop <- function(sub.id, which.DMCC, use.runs, do.shift) {     
  # sub.id <- "130114"; which.DMCC <- 2; use.runs <- c("Bas1", NA, "Pro1", "Pro2", "Rea1", "Rea2");
  error.str <- "";   # empty string for returning error messages
  
  # build a list of the input eprime files, into the all.ins object so they only need to be read in once.
  # same as in "D:\gitFiles_ccplabwustl\R01\Jo\for800msecTR\knitr\singleSubSummary\singleSubSummary.rnw"
  all.ins <- list(Bas1=NA, Bas2=NA, Pro1=NA, Pro2=NA, Rea1=NA, Rea2=NA);
  for (ssid in 1:length(session.ids)) {
    for (rid in 1:2) {     # ssid <- 3; rid <- 1;
      if (length(which(use.runs == paste0(sess.ids[ssid], rid))) == 1) {   # run should be good (listed in use.runs)
        if (wustl.box == TRUE) {
          fname <- paste0('"', sub.id, "_", session.ids[ssid], "_Stroop", sess.ids[ssid], "_run", rid, '.txt_raw"'); 
          tmp <- box_search(fname, type='file', file_extensions='csv', ancestor_folder_ids=folder.num); 
          if (length(tmp) == 1) { 
            in.tbl <- box_read_csv(tmp[[1]]$id);  
            in.tbl$Flicker.OnsetTime <- as.numeric(in.tbl$Flicker.OnsetTime);
            in.tbl$Fixation.OnsetTime <- as.numeric(in.tbl$Fixation.OnsetTime);
            in.tbl$FixationFinal.OnsetTime <- as.numeric(in.tbl$FixationFinal.OnsetTime);
            in.tbl$getready.OnsetTime <- as.numeric(in.tbl$getready.OnsetTime);
            in.tbl$Stimuli.OnsetTime <- as.numeric(in.tbl$Stimuli.OnsetTime);
            
            all.ins[[paste0(sess.ids[ssid], rid)]] <- in.tbl;
          } 
          if (length(tmp) == 0) { error.str <- c(error.str, paste("missing:", fname)); }
          if (length(tmp) > 1) { error.str <- c(error.str, paste("found > 1 file named", fname)); }
        } else {
          fname <- paste0(in.path, sub.id, "_", session.ids[ssid], "_Stroop", sess.ids[ssid], "_run", rid, ".txt_raw.csv"); 
          if (file.exists(fname)) { 
            in.tbl <- read.csv(fname); 
            all.ins[[paste0(sess.ids[ssid], rid)]] <- in.tbl; 
          } else { error.str <- c(error.str, paste("missing:", fname)); }
        }
      }
    }
  }
  
  # make the events files for all trials
  for (ssid in 1:length(session.ids)) {   # ssid <- 2;
    fout <- file(paste0(out.path, sub.id, "_Stroop_", session.ids[ssid], "_allTrials.txt"), 'wt');  # start an empty file for the onsets
    for (rid in 1:length(run.ids)) {    # rid <- 1;
      if (length(all.ins[[paste0(sess.ids[ssid], rid)]]) > 1) {   
        in.tbl <- all.ins[[paste0(sess.ids[ssid], rid)]];   # just to simplify the code
        start.value <- in.tbl$scanstart.RTTime[1];    # value which needs to be subtracted from all events onsets for the true onsets.
        if (is.na(start.value) | start.value < 1000) { stop("invalid start.value"); }
        inds <- which(!is.na(in.tbl$TrialType));
        onsets <- (in.tbl$Stimuli.OnsetTime[inds] - start.value)/1000;  # /1000 to convert to seconds
        if (do.shift) onsets <- onsets - (onsets %% 1.2);   # shift to closest previous TR.
        cat(paste0(paste(onsets, collapse=" "), " \n"), file=fout);     # add to the file with onsets for just this session
      } else {
        cat("*\n", file=fout);     # add to the file with onsets for just this session
      }
    }
    close(fout); unlink(fout);     # close the file we just wrote for this session and trial type
  }
  
  
  # make the events files for the Con and InCon TrialType trials, separately by LWPC: PC50 all sessions. bias = MC for bas; bias = MI for pro&rea; buff = MC for rea, cong.
  for (lbl in c("PC50", "bias", "buff")) {   # lbl <- "PC50";
    for (ttype in c("InCon", "Con")) {   # need one onset file for each $TrialType (Con and InCon)   # ttype <- "Con";
      for (ssid in 1:length(session.ids)) {   # ssid <- 2;
        if (session.ids[ssid] == "reactive") { need.length <- 120; } else { need.length <- 108; }   # ready for error-checking ... fixed number of trials.
        make.file <- TRUE;  # start off this loop
        if (session.ids[ssid] == "baseline" & lbl == "buff") { make.file <- FALSE; }
        if (session.ids[ssid] == "proactive" & lbl == "buff") { make.file <- FALSE; }
        if (session.ids[ssid] == "reactive" & lbl == "buff" & ttype == "InCon") { make.file <- FALSE; }
        if (make.file == TRUE) {  # this session has events of this type, so read in the (converted) eprime output and make an onset file
          fout <- file(paste0(out.path, sub.id, "_Stroop_", session.ids[ssid], "_", lbl, ttype, ".txt"), 'wt');  # start an empty file for the onsets
          for (rid in 1:length(run.ids)) {    # rid <- 1;
            if (length(all.ins[[paste0(sess.ids[ssid], rid)]]) > 1) {   
              in.tbl <- all.ins[[paste0(sess.ids[ssid], rid)]];   # just to simplify the code
              if (length(which(in.tbl$Procedure == "StroopTrialPROC")) != need.length) { stop(paste("not right number of trials:", fname)); }  # here's the trial error-checking
              start.value <- in.tbl$scanstart.RTTime[1];    # value which needs to be subtracted from all events onsets for the true onsets.
              if (is.na(start.value) | start.value < 1000) { stop("invalid start.value"); }
              if (lbl == "PC50") { inds <- which(in.tbl$TrialType == ttype & in.tbl$LWPC == "PC50"); }
              if (lbl == "bias" & session.ids[ssid] == "baseline") { inds <- which(in.tbl$TrialType == ttype & in.tbl$LWPC == "MC"); }
              if (lbl == "bias" & session.ids[ssid] != "baseline") { inds <- which(in.tbl$TrialType == ttype & in.tbl$LWPC == "MI"); }
              if (lbl == "buff" & session.ids[ssid] == "reactive") { inds <- which(in.tbl$TrialType == ttype & in.tbl$LWPC == "MC"); }
              if (length(inds) < 1) { stop(paste("no trials?", fname)); }
              onsets <- (in.tbl$Stimuli.OnsetTime[inds] - start.value)/1000;
              if (do.shift) onsets <- onsets - (onsets %% 1.2);   # shift to closest previous TR.
              cat(paste0(paste(onsets, collapse=" "), " \n"), file=fout);     # add to the file with onsets for just this session
              rm(inds);  # so can't accidently use them in the next loop.
            } else {
              cat("*\n", file=fout);     # add to the file with onsets for just this session
            }
          }
          close(fout); unlink(fout);     # close the file we just wrote for this session and trial type
        }
      }
    }
  }
  
  
  # block onset AND offset times; end the blocks at $Procedure == FixationGetReadyPROC, $getready.OnsetTime
  for (ssid in 1:length(session.ids)) {   # ssid <- 2;
    fout <- file(paste0(out.path, sub.id, "_Stroop_", session.ids[ssid], "_blockONandOFF.txt"), 'wt');  # start an empty file for the onsets
    for (rid in 1:length(run.ids)) {    # rid <- 1;
      if (length(all.ins[[paste0(sess.ids[ssid], rid)]]) > 1) {   
        in.tbl <- all.ins[[paste0(sess.ids[ssid], rid)]];   # just to simplify the code
        start.value <- in.tbl$scanstart.RTTime[1];
        inds <- which(in.tbl$Procedure == "FixationGetReadyPROC");   # same as above: block onsets
        if (length(inds) != 3) { stop(paste("not 3 FixationGetReadyPROC:", fname)); }
        onsets <- (in.tbl$getready.OnsetTime[inds] - start.value)/1000;  # Nick says use OnsetTime, not StartTime; /1000 to convert to seconds
        offset.1 <- (in.tbl$Fixation.OnsetTime[inds[2]] - start.value)/1000;
        offset.2 <- (in.tbl$Fixation.OnsetTime[inds[3]] - start.value)/1000;
        offset.3 <- (in.tbl$FixationFinal.OnsetTime[which(in.tbl$Procedure == "FixationOnlyPROC")] - start.value)/1000; # last coded differently
        onsets <- sort(c(onsets, offset.1, offset.2, offset.3));  # sort so onsets and offsets nicely arranged
        onsets <- onsets - (onsets %% 1.2);   # shift to closest previous TR.
        cat(paste0(paste(onsets, collapse=" "), " \n"), file=fout);     # add to the file with onsets for just this session
      } else { cat("*\n", file=fout); }    # add to the file
    }
    close(fout); unlink(fout);
  }
  
  
  # make the block timing files: start the blocks at onset of first trial in block; duration until start of last trial in block.
  for (ssid in 1:length(session.ids)) {   # ssid <- 1;
    fout <- file(paste0(out.path, sub.id, "_Stroop_", session.ids[ssid], "_block.txt"), 'wt'); 
    for (rid in 1:length(run.ids)) {    # rid <- 1;
      if (length(all.ins[[paste0(sess.ids[ssid], rid)]]) > 1) {   
        in.tbl <- all.ins[[paste0(sess.ids[ssid], rid)]];   # just to simplify the code
        start.value <- in.tbl$scanstart.RTTime[1];         # value which needs to be subtracted from all events onsets for the true onsets.
        if (is.na(start.value) | start.value < 1000) { stop("invalid start.value"); }          
        inds <- which(in.tbl$Procedure == "FixationGetReadyPROC");   # line before each block
        if (length(which(in.tbl$Procedure[inds+1] == "StroopTrialPROC")) != 3) { stop(paste("not expected row ordering", fname)); }
        if (in.tbl$Procedure[inds[2]-1] != "StroopTrialPROC" | in.tbl$Procedure[inds[3]-1] != "StroopTrialPROC") { stop(paste("not expected row ordering", fname)); }
        
        # find the first and last trial of each block, by looking relative to the FixationGetReadyPROC inds (which mark start of each block)
        onsets <- (in.tbl$Stimuli.OnsetTime[inds+1] - start.value)/1000;  # inds+1 to get first TRIAL of each block
        onsets <- onsets - (onsets %% 1.2);   # shift to closest previous TR.
        
        offset.1 <- (in.tbl$Stimuli.OnsetTime[inds[2]-1] - start.value)/1000;   # end of first block, in s
        offset.2 <- (in.tbl$Stimuli.OnsetTime[inds[3]-1] - start.value)/1000;
        offset.3 <- (in.tbl$Stimuli.OnsetTime[which(in.tbl$Procedure == "FixationOnlyPROC")-1] - start.value)/1000;  # last block, look for final fixation
        offset.1 <- offset.1 - (offset.1 %% 1.2);   # shift to closest previous TR.
        offset.2 <- offset.2 - (offset.2 %% 1.2); 
        offset.3 <- offset.3 - (offset.3 %% 1.2); 
        if (onsets[2]-offset.1 < 30 | onsets[3]-offset.2 < 30 | onsets[1] < 30) { stop(paste("breaks too short:", fname)); }
        # duration of each block is offset - onset
        cat(paste0(onsets[1], ":", (offset.1-onsets[1]), " ", onsets[2], ":", (offset.2-onsets[2]), " ", onsets[3], ":", (offset.3-onsets[3]), "\n"), file=fout);   
      } else { 
        cat("*\n", file=fout);     # * if run missing
      } 
    }
    close(fout); unlink(fout);     # close the file we just wrote for this session and type
  }
  
  return(list(err.str=error.str, eprimes=all.ins));
}




wustl.box <- FALSE;   # if a DMCC team member is running this code (downloads the eprime csv input from box)

sub.id <- "132017";  #  subject id. (in quotes, even if the subject id is all numbers.)
which.DMCC <- 2;    # which.DMCC <- 3;  # DMCC phase number (e.g. which.DMCC <- 3 for the second wave of scans, DMCC_Phase3)

setwd("C:/Users/mcf/Box/global/proj/ub55")
out.path <-  file.path(getwd(), "out", "glms", sub.id, "INPUT_DATA", "Stroop");     # path to a local directory where the files will be written (in a sub.id subdirectory).

# download evtFileMaker.R from https://github.com/ccplabwustl/dualmechanisms/tree/master/preparationsAndConversions/eprime/
code.fname <- "d:/gitFiles_ccplabwustl/dualmechanisms/preparationsAndConversions/eprime/evtFileMaker.R";  # local path to evtFileMaker.R

# in.path isn't used if wustl.box <- TRUE;  
if (wustl.box == FALSE) { in.path <- "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Raw_Data/132017"; }  # path to eprime txt_raw.csv files for this person (written by TEMPLATE_convertEprime.R)



do.Stroop(sub.id, which.DMCC, c("Bas1", "Bas2"))
