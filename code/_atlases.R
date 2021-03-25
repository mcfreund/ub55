
## atlases ----

hcp <- list(
  L  = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii")
  ),
  R = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii")
  )
)


if (nodename == "ccp-freund") {
  
  hcp.vinf <- list(
    L  = readGIfTI(
      file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii")
    ),
    R = readGIfTI(
      file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii")
    )
  )
  
}



parcellation <- mikeutils::read_atlas("schaefer400", path.atlas = dir.schaefer)

if (nodename == "ccplinux1") {
  
  schaefer10k <-
    c(
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
    )

}


over <- list(
  L = parcellation$atlas[1:(nrow(parcellation$atlas) / 2)], 
  R = parcellation$atlas[(nrow(parcellation$atlas) / 2):nrow(parcellation$atlas)]
)

schaefer <- list(
  L = read_gifti2matrix(
    file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_L.label.gii")
  ) %>% c,
  R = read_gifti2matrix(
    file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_R.label.gii")
  ) %>% c
)


