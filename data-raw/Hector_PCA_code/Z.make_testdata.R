#### Make data for PCA tests

datadir <- 'PCA-data'

test_pca <- readRDS(file.path(datadir, 'pc-emissCC-esmHistorical-esmrcp85.rds'))
tc <- dplyr::bind_rows(lapply(file.path(datadir,
                                        c('emissCC-esmHistorical.rds','emissCC-esmrcp85.rds')),
                              readRDS))

## trim the test data to just 20 PCs
ncomp <- 20
test_pca$sdev <- test_pca$sdev[1:ncomp]
test_pca$rotation <- test_pca$rotation[ , 1:ncomp]

p1 <- rep(0,ncomp) ; p1[1] <- 1
p2 <- rep(0,ncomp) ; p2[1] <- 1; p2[2] <- 0.5
plist <- list(p1,p2)
test_runid <- min(tc$runid[tc$runid > length(plist)])
tc <- dplyr::filter(tc, runid==test_runid)

test_climdata <- lapply(seq_along(plist),
                        function(i) {
                            p <- plist[[i]]
                            cd <- hectorcal::reconstruct_climate(p, test_pca)
                            cd$runid <- i
                            cd
                        })
test_climdata <- dplyr::bind_rows(test_climdata, tc)

saveRDS(test_pca, file.path(datadir, 'test_pca.rds'))
saveRDS(test_climdata, file.path(datadir, 'test_climdata.rds'))


