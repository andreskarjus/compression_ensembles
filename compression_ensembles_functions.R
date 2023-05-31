############################### #
# Functions to accompany the script file for the paper 
# "Compression ensembles quantify aesthetic complexity and the evolution of visual art"
# by Karjus et al 2023
# This software is provided as-is, without any implication of warranty.
# If you make use of this software or the precomputed datasets, kindly cite the paper.
#
# This is the functions file, intended to be sourced from the script file.
############################## #


p=c("MASS", "magick", "caret", "data.table", "mgcv", "Hmisc", "tidyverse", "entropy", "gifski",  "reshape2", "parallel", "Metrics", "patchwork", "colorspace", "scales", "umap") 
# mass needs to precede dplyr/tidyverse, so dplyr::select can take precedence
install.packages(setdiff(p, rownames(installed.packages())))  # installs if missing
print(sapply(p, require, character.only = TRUE))



do_pca = function(comps, dat2, dontsign=NULL, trim=T, exclude=NULL, corvar="png_x1"){
  compvecs = comps %>% .[sapply(.,length)>1]  # unpack chunks, remove fails
  print(paste(length(compvecs), "vectors out of", length(comps)))
  fps = sapply(compvecs, function(y) attr(y, "fp"))
  compvecs = do.call(rbind, compvecs) 
  if(!is.null(exclude)){
    compvecs = compvecs[, !(colnames(compvecs) %in% exclude) ]
  }
  compvecs = compvecs[,  !duplicated(colnames(compvecs))]
  # summary(pr)$importance[3,] %>% {plot(.); abline(v=which(.==1)[1])}
  var0 = which(apply(compvecs, 2, var)==0)
  if(length(var0)>0){
    print(paste("Zero-variance columns detected, removing:", paste(names(var0), collapse=" ") ))
    compvecs = compvecs[,-var0]
  }
  
  pr0 = prcomp(compvecs, scale. = T, center = T)
  
  loads = pr0$rotation %>% abs()# %>% apply(.,2, function(x) x/sum(x) ) # in global for debug
  #vrs = t(t(loads)/rowSums(t(loads)))*100  # relative contibution of each var
  
  topvars=list()
  for(i in 1:ncol(loads)){
    v = sort(loads[,i], decreasing = T)
    topvars[[i]]=(head(v,5))
    # v2 = which.max(loads[,i])
    # cr = cor(pr0$x[,i], compvecs[, v2])
    # if(cr < 0){
    #   pr0$x[,i] = -1 * pr0$x[,i] # change sign if needed
    # } 
    
    # better: correlate with baseline gif, and make the sign of the PC the same
    # this way most PCs ordered as simple<-->complex
    cr = cor(pr0$x[,i], compvecs[, corvar])
    if(cr < 0){
      pr0$x[,i] = -1 * pr0$x[,i] 
    }
  }
  
  # debug
  # cr=c()
  # for(i in 1:ncol(pr0$x)){
  #   cr[i] = abs(cor(compvecs[,"gif"], pr0$x[,i]))
  # }
  # cr=c()
  # for(i in 1:ncol(compvecs)){
  #   cr[i] = abs(cor(compvecs[,i], pr0$x[,1]))
  # }
  # colnames(compvecs)[which.max(cr)]
  # max(cr)
  
  
  pr0 <<- pr0  # for plots and debug
  loads <<- loads
  #vrs <<-vrs
  compvecs <<- compvecs
  topvars <<- topvars
  
  if(trim){
    pr = pr0$x[,1:which(summary(pr0)$importance[3,]==1)[1] ]
  } else {
    pr = pr0$x
  }
  
  dat = dat2[dat2$Path2 %in% fps, ] %>% ungroup()
  if(!all(fps == dat$Path2)){
    stop("Metadata doesn't match vectors, something is off")
  } # double check; with the test set of multiple langs it reuses vector automatically
  dat = cbind(dat, as.data.frame(pr))
  print(paste(nrow(dat), "objects ready"))
  #return(list(dat=dat, loads=loads ))
  return(dat)
}





dataloader=function(cfiles, datpath, corvar="compress_gif_x1", hic=F){
  load(datpath)
  if(!hic){
    dat = dat %>% select(author_name, painting_name,Date, Date2,  
                         Genre, Style, "Art Movement", Field, Media, 
                         Dimensions, Path2, Pathb, image_url) %>% ungroup()
  }
  
  compvecs=list()
  for(i in cfiles){
    load(i)  # comps object
    compvecs=c(compvecs, 
               comps %>% 
                 unlist(recursive = F) %>% 
                 .[sapply(.,length)>1])
    rm(comps)
  }
  
  fps = sapply(compvecs, function(y) attr(y, "fp")) %>% basename()
  
  # normalize later; also don't need 15th digit precision?
  compvecs = do.call(rbind, compvecs) # %>% round(10)  
  compvecs = compvecs[,  !duplicated(colnames(compvecs))]
  
  dat = dat[dat$Pathb %in% fps,,drop=F ] 
  dat = dat[match(fps, dat$Pathb),,drop=F ] 
  if(!all(fps == dat$Pathb) | nrow(compvecs) != nrow(dat) | any(duplicated(fps)) ){
    stop("Metadata doesn't match vectors, something is off")
  } # double check; with the test set of multiple langs it reuses vector automatically
  
  # summary(pr)$importance[3,] %>% {plot(.); abline(v=which(.==1)[1])}
  var0 = which(apply(compvecs, 2, var) %>%  {is.na(.)|.==0} )
  if(length(var0)>0){
    print(paste("Zero-variance columns detected, removing:", paste(names(var0), collapse=" ") ))
    compvecs = compvecs[,-var0]
  }
  
  pr0 = prcomp(compvecs, scale. = T, center = T) # 
  loads = pr0$rotation %>% abs()
  topvars=list()
  for(i in 1:ncol(loads)){
    v = sort(loads[,i], decreasing = T)
    topvars[[i]]=names(head(v,5))
    # correlate with baseline gif, and make the sign of the PC the same
    # this way most PCs ordered as simple<-->complex
    cr = cor(pr0$x[,i], compvecs[, corvar])
    if(cr < 0){
      pr0$x[,i] = -1 * pr0$x[,i] 
    }
  }
  
  
  colnames(compvecs) = paste0("TR_", colnames(compvecs))
  dat = cbind(dat, as.data.frame(compvecs), as.data.frame(pr0$x %>% round(10)))
  attr(dat, "loads")=pr0$rotation %>% abs()
  attr(dat, "topvars") = topvars
  
  # data fixes
  if(!hic && "author_name" %in% colnames(dat) && "Georgia O Keeffe" %in% dat$author_name){
    dat$author_name[dat$author_name=="Georgia O Keeffe"]="Georgia O'Keeffe" 
  }
  
  
  print(paste(nrow(dat), "objects ready"))
  return(dat)
  
}




dataloader_evals = function(compvecs, dat){
  fps = sapply(compvecs, function(y) attr(y, "fp")) %>% basename()
  compvecs = do.call(rbind, compvecs) %>% scale()  # scale here, but won't do pca for this
  compvecs = compvecs[,  !duplicated(colnames(compvecs)) ]
  
  dat = dat[dat$Pathb %in% fps, ] 
  
  if(!all(fps == dat$Pathb)  | any(duplicated(fps)) ){ # test reuses, ok
    stop("Metadata doesn't match vectors, something is off")
  } # double check; with the test set of multiple langs it reuses vector automatically
  
  # summary(pr)$importance[3,] %>% {plot(.); abline(v=which(.==1)[1])}
  var0 = which(apply(compvecs, 2, var) %>%  {is.na(.)|.==0} | apply(compvecs, 2, function(x)length(unique(x))) < 10   )  # extra filter as multipic otherwise gets a few variables with 1 unique value (due to something in how the images are produced I guess), and fractals has 3 not 1 contrast sd values; but these would throw off the classifiers.
  
  if(length(var0)>0){
    print(paste("Low or zero-variance columns detected, removing:", paste(names(var0), collapse=" ") ))
    compvecs = compvecs[,-var0]
  }
  
  colnames(compvecs) = paste0("PC_", colnames(compvecs))
  return(cbind(dat, as.data.frame(compvecs)) )  # reuses df multiple times to match
}


do_all_evals = function(projectfolder, multipicdata, multipicnorms, fractaldata, fractalnorms){
  print("doing multipic")
  params = list(
    nfree = 0
    ,projectfolder = projectfolder
    ,corpusfolder = file.path(projectfolder, multipicdata) 
    ,datfolder = file.path(projectfolder, multipicnorms) 
    #"C:/Users/Andres/korpused/artcorpora/norms/PicPsy/pics"
    ,pixeltotal=300*300   # actual size of multipic ones; 1 is smaller
    ,pixelsize=300          
    ,coefs = c(1,0.4)
    ,extracoefs=c( 0.2, 0.1)
    ,minsize=0.5
  )
  
  f = list.files(params$datfolder, full.names = T, pattern = "csv$")
  dat2 = lapply(f, function(ff)
    read.csv(file = ff, sep=";",dec = ",", encoding = "UTF-8" )  %>% select(PICTURE, VISUAL_COMPLEXITY) %>% filter(!is.na(VISUAL_COMPLEXITY)) %>%  mutate(language=gsub("_MultiPic_CSV.csv$","", basename(ff))) %>% mutate(Path2 = file.path(params$corpusfolder, paste0(PICTURE, ".png")) )
  )  %>% do.call(rbind, .) %>% mutate(Pathb = basename(Path2))
  # same set of images for all languages, so only need to calculate once
  comps = docompressions_multicore(dat2 %>% filter(language=="English") %>% select(Path2), params, pathvar="Path2")
  print(paste(Sys.time(), "compressions done"))
  # dat = do_pca(comps, dat2)
  dat = dataloader_evals(comps, dat2 ) # error sapply(compvecs, function(y) attr(y, "fp")) %>% basename()
  
  datl=dat %>% mutate(VISUAL_COMPLEXITY = rescale(VISUAL_COMPLEXITY, to=c(1,5))) %>% 
    mutate(language=paste0("Multipic_", language)) %>% mutate(gif=PC_compress_gif_x1) 
  
  save(datl,comps,  file=file.path(params$projectfolder, "normtest2_gif.RData"))
  
  print("doing fractals")
  params = list(
    nfree = 0
    ,projectfolder = projectfolder
    ,corpusfolder = file.path(projectfolder, fractaldata) 
    ,datfolder = file.path(projectfolder, fractalnorms) 
    #"C:/Users/Andres/korpused/artcorpora/norms/PicPsy/pics"
    ,pixeltotal=380*380   
    ,pixelsize=380          
    ,coefs = c(1, 0.4)
    ,extracoefs=c(0.2,0.1)
    ,minsize=0.5
  )
  f = list.files(params$datfolder, full.names = T, pattern = "csv$")
  dat2 = read.csv(file = file.path(params$datfolder, "fractals_norms.csv"), sep=",",dec = ".", encoding = "UTF-8" )  %>% select(fractal,zip.ratio, complexity) %>% filter(!is.na(complexity)) %>% mutate(Path2 = file.path(params$corpusfolder, paste0(fractal, ".bmp"))) %>% mutate(Pathb = basename(Path2))
  comps = docompressions_multicore(dat2 %>% select(Path2), params, pathvar="Path2")
  #dat = do_pca(comps, dat2)
  dat = dataloader_evals(comps, dat2)
  datf = dat %>% rename(VISUAL_COMPLEXITY=complexity) %>% mutate(gif=PC_compress_gif_x1)  %>% 
    mutate(VISUAL_COMPLEXITY = rescale(VISUAL_COMPLEXITY, to=c(1,5)))
  save(datf, comps,  file=file.path(params$projectfolder, "normtest_fractals2.RData"))
  
  
  #### load test data, test and and plot ####
  print("doing tests")
  #compvecs[, grep("png|gif|jpeg|tiff|svg", colnames(compvecs))] %>% cor() %>%  reshape2::melt() %>% ggplot(aes(Var1, Var2, fill=value))+geom_tile()+scale_fill_viridis_c()
  
  # crossvalidate for all languages, report error, also for fractals, compare gif vs pca vs artmodel pca vs lite model with top variables only
  #topvars[c(1,9,8,2,3)] %>% sapply(function(x) paste(x[1:2], collapse="  ") ) 
  
  langs = c(unique(datl$language), "Fractals")
  nl = length(langs)
  # rmse=data.frame(lang=rep(langs, each=25), r=NA)
  rmse=vector("list", nl); names(rmse)=langs
  
  r2=data.frame(lang=rep(langs,each=3), r=NA, rtype=rep(c("top5", "full", "gif"), nl))
  ii=1; ii2=1
  for(l in  1:nl){
    if(langs[l] == "Fractals"){
      d = datf
    } else {
      d = datl %>% 
        filter(language==langs[l])
    }
    
    d = d %>% select( -which(grepl("PC", colnames(.)) & 
                               apply(d,2,function(x) length(unique(x)))<10 ) )
    
    rs=c()
    n=d %>% select(starts_with("PC")) %>% ncol()
   
    ii2=ii2+1
    m = (lm(VISUAL_COMPLEXITY~., 
            data=d %>% select(starts_with("PC") | starts_with("VISUAL") )))
    r2[ii2,"r"] = summary(m)$adj.r.squared
    
   
    
    ii2=ii2+1
    m = summary(lm(VISUAL_COMPLEXITY~gif, data=d))
    r2[ii2,"r"] = m$adj.r.squared
    ii2=ii2+1
    
    for(i in 1:50){
      dd=d %>% select(VISUAL_COMPLEXITY, starts_with("PC") )
      #dd = d %>% select_if(colnames(.) %in% c("VISUAL_COMPLEXITY", names(sort(rs, decreasing = T))[1:10] ))
      n = nrow(dd)
      splits=round(n*0.9)
      tr = sample(1:n, splits)
      te = setdiff(1:n, tr)
      #m = (lm(VISUAL_COMPLEXITY~PC1+PC2+PC3+PC4+PC15, data=dd[tr,]))
      m = lm(VISUAL_COMPLEXITY~., data=dd[tr,])
     
      rmse[[l]] = c(rmse[[l]], unname( abs(predict(m, dd[te,]) - dd[te,"VISUAL_COMPLEXITY"])) )
      
      # ii=ii+1
    }
    
  }
  save(datl,datf,rmse,r2, file=file.path(params$projectfolder, "human_evals.RData"))
  
  sapply(rmse, median)  %>% print()
  r2 %>% group_by(lang, rtype) %>% summarize(mean(r)) %>% print()
  #datl$VISUAL_COMPLEXITY %>% sd
  #rmse<<-rmse
  r2<<-r2
  return(rmse)
}


# simple shuffle-one permute won't work if correlated features; will average binary classifiers instead, using z-values as a rough estimate of variable importance; can sort the matrix later too if needed for a clearer visualization (the number of vars really matters). 
# also can easily get rid of some correlated vars by removing the small-transformer results.
approx_importance = function(datx, nclass, ntest, nruns=100, fullcomb=F, continuous=F){
  #apply(datx[,-1], 2, sd) %>% sort() %>% head
  # remove couple of near-constant variables, would break the lda later:
  datx = datx[,c("xvar", 
                 apply(datx[,-1],2, function(x) length(unique(x))) %>% 
                   .[.>(nrow(datx)/10)] %>% names() )]
  cl=colnames(datx)[-1]
  if(!continuous){
    lv = unique(datx$xvar)
    datx = datx %>% mutate(xvar=as.factor(xvar)) #%>% as.data.frame()
    datx[,-1] = scale( datx[,-1]) # normalize (re-enabled as input is raw)
    if(fullcomb){
      # for many-multiclass, try all combos (or half at least)
      com = combn(lv,2) %>% .[,sample(1:ncol(.))] %>% .[,1:(ncol(.)/2)]
      nruns=ncol(com)
    }
  }
  
  res = matrix(NA, ncol = nruns, nrow=ncol(datx)-1, # first is class column
               dimnames=list(cl,list()))
  for(j in 1:nruns){
    if(continuous){
      x = datx  %>% 
        ungroup() %>% 
        slice(sample(1:n(), nclass )) 
      res[,j]= lm(xvar~., x) %>% varImp() %>% .[,1]
    } else {
      if(fullcomb){
        lvx=com[,j,drop=T]
      } else {
        lvx = sample(lv, 2)
      }
      x = datx  %>% 
        ungroup() %>% 
        filter(xvar %in% lvx) %>%  # build binary classifier for rndm pair
        droplevels() %>% 
        group_by(xvar) %>%  
        slice(sample(1:n(), pmin(nclass, n()) )) # %>% # use some random n subset of data
      #mutate(test = row_number() %in% sample(1:n(), ntest))  # subset of that for test
      #xtest =  x %>% filter(test) %>% select(-c(test))# %>% as.data.frame()
      # res[j,i] = 
      res[,j]= glm(xvar~., x, family="binomial") %>% varImp() %>% .[,1]
      # %>% arrange(Overall)
      # collect those, take mean, done.
    }
  }
  res = rowMeans(res, na.rm=T) %>% sort(decreasing = T) %>% names()
  if(!continuous){
    x3 = head(res, 3)
    res2 = matrix(NA, ncol=20, nrow=3, dimnames=list(x3, list()))
    # the top ones: compare which improves model most, order by that (affects the plot looks the most so should be as accurate as possible; but doing it for all var combos would take ages)
    for(v in seq_along(x3)){
      for(j in 1:20){
        x = datx[,c("xvar", x3[v])]  %>% 
          group_by(xvar) %>%  
          slice(sample(1:n(), pmin(nclass, n()) )) %>% # use some subset of data
          mutate(test = row_number() %in% sample(1:n(), ntest))  # subset of that for test
        xtest =  x %>% filter(test) %>% select(-c(test))# %>% as.data.frame()
        res2[v,j] = (predict(lda(xvar~., 
                                 x %>% filter(!test) %>% select(-c(test)), method="mle"),
                             xtest)$class == xtest$xvar) %>% mean(na.rm=T)
      }
    }
    res2 = rowMeans(res2, na.rm=T) %>% sort(decreasing = T) %>% names()
    res2 = c(res2[1], setdiff(res, res2[1]))
  } else {
    res2=res
  }
  
  # move collinear variables to the back
  cr = abs(cor(datx[,res2], datx[,res2]))
  diag(cr)=0
  ok=setNames(rep(0L,length(res2)), res2)
  for(i in 1:(length(ok)-1) ){
    if(ok[i]==0){
      cr2 = as.integer(cr[i,]>0.6)
      cr2 = cr2 + (as.integer(cr[i,]>0.8) * 100)
      cr2[1:i] = 0
      ok = ok + cr2
    }
  }
  ok=ok[ok<100]
  return(names(sort(ok))) # preserves order where values equal
}



detection_accuracy = function(datx,xvars, 
                              nclass=2500, ntest=250,maxtrain=2000, nruns=100, 
                              trainseq=c(10,100,1000),
                              kappa=F, debug=F, prefix=NULL, doconfusion=T,
                              baselinevar="TR_compress_gif_x1",
                              allvarcombos=T){
  # input data should have min ntest+maxtrain*n groups
  
  # reorder to show gif baseline in plot; doesn't affect full model results
  if(!is.null(baselinevar)){
    xvars = c(baselinevar, setdiff(xvars, baselinevar))
  }
  
  datx = datx[, c("xvar", xvars)] # this now has constant and collinear ones removed
  if(allvarcombos){
    varseq = 1:length(xvars)
    vn = xvars[varseq] %>% gsub("^TR_|_x1$", "",.)
  } else {
    varseq = c(1:10, 30, 50, length(xvars))
    vn = xvars[varseq] %>% gsub("^TR_|_x1$", "",.) # semimanual hacky labels 
    vn[2:length(vn)] = paste0("+",vn[2:length(vn)])
    vn[length(vn)-2] = "+next 20 transforms"
    vn[length(vn)-1] = "+another 20 transforms"
    vn[length(vn)] = paste0("+all ", length(xvars), " transforms")
  }
  
  
  if(debug){
    varseq=length(xvars)
    vn=paste0("+all ", length(xvars), " transforms")
  }
  
  datx = datx %>% mutate(xvar=as.factor(xvar)) #%>% as.data.frame()
  # for kappa:
  r = 1/length(unique(datx$xvar)) # assuming roughly equal size classes
  r1=1-r
  lastkappa=NULL # if not kappa for all then add as attribute for last
  datx[,-1] = scale( datx[,-1]) # normalize as input is now raw transforms  
  res = matrix(NA, 
               nrow=length(varseq), # first is class column
               ncol = length(trainseq), 
               dimnames=list(vn, as.character(trainseq)))
  # for confusion matrix (if needed):
  cm = table(levels(datx$xvar),levels(datx$xvar)) %>% {.[]=0;.} 
  
  library(parallel)
  cl = makeCluster(detectCores())
  tryCatch({
    clusterEvalQ(cl, suppressPackageStartupMessages(
      c(library("MASS"), library("dplyr")))) # dplyr last bc of select
    clusterExport(cl, c("nclass", "ntest", "maxtrain", "xvars", "varseq","trainseq", "datx", "doconfusion"), envir=environment() )
    
    
    for(v in seq_along(varseq)){
      for(ne in seq_along(trainseq)){
        # tmp=vector("list", nruns )
        # for(j in 1:nruns){
        clusterExport(cl, c( "v", "ne"), envir=environment() )
        tmp = parLapply(cl, 1:nruns, function(j) {
          
          x = datx[,c("xvar", xvars[1:varseq[v]]) ] %>% 
            group_by(xvar) %>%  
            slice(sample(1:n(), pmin(nclass, n()) )) %>% # use some subset of data
            mutate(test = row_number() %in% sample(1:n(), ntest))# %>%   # subset of that for test
          # if(v==length(varseq) & ne==length(trainseq)){
          #   x = x %>% {data.frame(xvar=.$xvar,test=.$test, prcomp(.[-1], scale. = T, center = T)$x )} %>% group_by(xvar)
          # }
          # 
          xtrain = x %>% filter(!test) %>% select(-c(test)) %>% slice(sample(1:n(), pmin(n(), trainseq[ne])))
          xtest =  x %>% filter(test) %>% select(-c(test))
          
          try({
            pr = predict(lda(xvar~., xtrain, method="mle"), xtest)$class
            # if(doconfusion){
            #   if(ne==length(trainseq) & v==length(varseq)){
            #     cm = cm + table(xtest$xvar, pr) # rows = truth, col=prediction
            #   }
            # }
          })
          #return((pr == xtest$xvar))
          return(list(pr=pr, xtest=xtest$xvar))
        }
        ) # end parLapply
        # tmp = unlist(tmp, F,F)
        gc(verbose=F)
        if(doconfusion){
          if(ne==length(trainseq) & v==length(varseq)){
            for(i in seq_along(tmp)){
              cm = cm + table(tmp[[i]]$xtest, tmp[[i]]$pr)
            }
            # cm = cm + table(xtest$xvar, pr) # rows = truth, col=prediction
          }
        }
        
        tmp = unlist(lapply(tmp, function(x) x$pr==x$xtest), F,F)
        
        if(kappa ){ # & !is.na(res[v,ne])
          # use the random success baseline as the prior expected prob for more accurate assessment
          res[v,ne] = ifelse(binom.test(c(table(tmp)), p = r)$p.value<0.05, 
                             mean(tmp, na.rm=T), NA)
          res[v,ne] = (res[v,ne]-r)/r1
        } else {
          res[v,ne] = mean(tmp, na.rm=T) # mean of replicates
          
          
        }
        if(ne==length(trainseq) & v==length(varseq)){
          se = sd(tmp, na.rm = T)/sqrt(length(tmp))
          if(!kappa){ # if rest is acc then add kappa for last one
            lastkappa = (res[v,ne]-r)/r1
          }
        }
        
      } # end loops
    }
    
  }, # end trycatch
  error=function(e){print(e)},
  warning=function(e){print(e)},  
  finally = stopCluster(cl) 
  )
  
  res=reshape2::melt(res, value.name = "acc") %>%
    rename(v=Var1, ne=Var2) %>% 
    mutate(ne=as.factor(ne)  )   %>% 
    group_by(ne) %>% mutate(nvars=row_number()) %>% ungroup()
  # v=factor(v,levels =  vn) 
  #levels(res$v) = paste0( rev(letters)[ 1:length(levels(res$v))] ,prefix, "_", levels(res$v) )  # was useful for matrix plotting only
  
  if(doconfusion){
    attr(res, "cm")=cm
  }
  attr(res, "se")=se
  attr(res, "kappa")=lastkappa
  return(res )
}


do_common_pca = function(dat, hdat){
  d1 = dat %>% select(starts_with("TR_", ignore.case = F))
  d2 = hdat %>% select(starts_with("TR_", ignore.case = F))
  d = rbind(d1[, intersect(colnames(d1), colnames(d2))],
            d2[, intersect(colnames(d1), colnames(d2))] )
  
  pr0 = prcomp(d, scale. = T, center = T) # 
  loads = pr0$rotation # %>% abs()
  topvars=list()
  for(i in 1:ncol(loads)){
    v = loads[,i][order( abs(loads[,i]), decreasing = T)]
    topvars[[i]]=(head(v,5))
    # correlate with baseline gif, and make the sign of the PC the same
    # this way most PCs ordered as simple<-->complex
    cr = cor(pr0$x[,i], d[, "TR_compress_gif_x1"])
    if(cr < 0){
      pr0$x[,i] = -1 * pr0$x[,i] 
    }
  }
  p = pr0$x
  l = list(
    art = cbind(dat[,1:8],Path2 = dat$Path2, p[1:nrow(dat),]),
    hic = cbind(hdat[,1:10],Date2=hdat$Date2,Path2 = hdat$Path2, p[(nrow(dat)+1):nrow(p),])
  )
  attr(l, "topvars")=topvars
  return(l)
}


do_mean_sd_space = function(dd, miny=1800,maxy=1990,minwin=10,maxwin=50,target=500, maxpc=3, nbin=NULL){
  dd = arrange(dd, Date2)
  # bin_equal = function(x, nbin = nbin) {
  #   breaks = quantile(x, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE)
  #   return(findInterval(x, breaks, all.inside = TRUE))
  # }
  ysize = rep(NA, 2000)
  xt = table(dd$Date2)[as.character(min(dd$Date2):max((dd$Date2)))]  %>% 
    {names(.)=min(dd$Date2):max((dd$Date2));.} %>% replace_na(0)
  win = frollsum(xt, c(minwin:maxwin), align = "center", fill = 0)  %>% do.call(cbind,.) %>% {.>=target} %>% cbind(T) %>% 
    {c(minwin:maxwin, maxwin)[apply(.,1,function(x) which(x)[1] )] } %>% { cbind(w= ., xt)}
  
  # dat %>% mutate(Move2 = `Art Movement`  %>% gsub("^([^,]+),.*", "\\1",.)) %>% group_by(Move2) %>% mutate(nm=n()) %>% filter(nm>200) %>% mutate(m=mean(Date2)) %>% ggplot(aes(y=reorder(Move2,m), x=Date2))+geom_violin()+geom_boxplot(width=0.2, outlier.colour = NA)
  
  # hicetnunc: will use days, counting from day 1
  dt = list()
  for(i in paste0("PC", 1:maxpc)){
    ysd=NA;ym=NA
    # ysd = sapply(seq(miny, maxy, 1), function(x) (dd %>% filter(Date2 %in% ( (x-win) :(x+win-1)))) %>% 
    #                .[,i,drop=T] %>% bootsd() )  %>% t()
    # colnames(ysd) = paste0("s", colnames(ysd))
    # ym = sapply(seq(miny, maxy, 1), function(x) (dd %>% filter(Date2 %in% ( (x-win) :(x+win-1)))) %>% 
    #               .[,i,drop=T] %>% bootmean() ) %>% t()
    # colnames(ym) = paste0("m", colnames(ym))
    
    # ymed = sapply(seq(miny, maxy, 1), function(x) (dd %>% filter(Date2 %in% ( (x-win) :(x+win-1)))) %>% 
    #                .[,i,drop=T] %>% quantile( probs=c(0.25,0.5,0.75)) )  %>% t() %>% as.data.frame()
    # colnames(ymed)=c("q1", "q2", "q3")
    
    # dt[[i]] = dd %>% mutate(bin = bin_equal(Date2, nbin)) %>% group_by(bin) %>% summarise(Date2=median(Date2), n=n(), value = quantile(!!as.symbol(i), c(0.025, 0.25,0.5,0.75, 0.975)), q= c(0.025, 0.25,0.5,0.75, 0.975) ) %>% ungroup() %>% pivot_wider(c(bin,n, Date2, value), names_from = q, names_prefix = "q") %>% mutate(pc=i)
    
    ydat = win[as.character(seq(miny, maxy, 1)), ]
    ytmp=tibble()
    for(y in seq(miny, maxy, 1) ){
      w = win[as.character(y), "w"]
      ymed = dd %>% filter(Date2 %in% ( (y-w) :(y+w)) ) %>% 
        .[,i,drop=T] %>% quantile( probs= c(0.025, 0.25,0.5,0.75, 0.975) )   %>%
        t() %>% as.data.frame()
      ytmp = rbind(ytmp, ymed)
    }
    colnames(ytmp)= paste0("q", c(0.025, 0.25,0.5,0.75, 0.975))
    
    # ymed = sapply(seq(miny, maxy, 1), function(x) (dd %>%
    #       filter(Date2 %in% ( (x-win) :(x+win-1))
    #              ) ) %>%
    # .[,i,drop=T] %>% quantile( probs= c(0.025, 0.25,0.5,0.75, 0.975) ) )  %>%
    #   t() %>% as.data.frame()
    #  colnames(ymed)= paste0("q", c(0.025, 0.25,0.5,0.75, 0.975))
    dt[[i]] =  cbind(ydat, ytmp, data.frame(Date2=seq(miny, maxy, 1), pc=i))
    
  }
  # y = sapply(seq(miny, maxy, 1), function(x) (dd %>% filter(Date2 %in% ( (x-win) :(x+win-1)))) %>% nrow())
  #dt[["x"]] = data.frame(year=seq(1800, 1990, 1),pc="x", sd=y/1000, mean=NA )
  dt = do.call(rbind, dt) # %>% filter(Date2 >= miny, Date2 < maxy)
  
  return(dt)
}


medplot = function(dots, dtart, r=c(-16,9), sub=""){
  update_geom_defaults("line", list(lineend = "round", linejoin="round"))
  notch1 = min(dtart$Date2) + ((max(dtart$Date2)-min(dtart$Date2))*0.05)
  notch2 = max(dtart$Date2) - ((max(dtart$Date2)-min(dtart$Date2))*0.05)
  g = ggplot( data=dtart, mapping=aes(x=Date2))+
    geom_hline(yintercept = 0, color="gray80")+
    rasterise(
      geom_point( mapping=aes(y=value), data=dots, alpha=1, size=0.01, shape=15, color=dots$col,fill=NA, position=position_jitter(height = 0, width=0.49)),
      dpi = 150, dev = "ragg_png"
    )+
    geom_line(aes(x=Date2,y=q0.025, color=w,alpha=w), color="gray70", lineend = "round", linejoin="round")+ #alpha=0.2)+
    geom_line(aes(x=Date2,y=q0.975, color=w,alpha=w), color="gray70", lineend = "round", linejoin="round")+ #alpha=0.2)+
    geom_line(aes(x=Date2,y=q0.25, color=w,alpha=w),  color="gray30", lineend = "round", linejoin="round")+ #alpha=0.5)+
    geom_line(aes(x=Date2,y=q0.75, color=w,alpha=w),  color="gray30", lineend = "round", linejoin="round")+ #alpha=0.5)+
    geom_line(aes(x=Date2,y=q0.5, color=w, alpha=w), color="black", lineend = "round", linejoin="round")+
    
    annotate("segment", x=-Inf, xend=notch1,y=c(-4.95,9.95),yend=c(-4.95,9.95), color="black", size=0.8)+ # compare-plots-line
    annotate("segment", x=notch2, xend=Inf,y=c(-4.95,9.95),  yend=c(-4.95,9.95), color="black", size=0.8)+
    
    #scale_color_gradientn(colors=c("black", "gray80"))+
    scale_alpha(range=c(0.05,1) %>% rev, limits=c(10,25))+
    #scale_fill_gradient(low="gray99", high="gray8")+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          #axis.text.y = element_blank(), #  element_text(size=7),
          axis.text.x = element_text(size=rel(1)),
          axis.text.y = element_text(size=rel(1)),
          axis.title=element_blank(),
          #axis.title.y = element_text(size=13, margin = margin(0,0,0,0)),
          axis.ticks.length = unit(1.5, "pt"),
          #axis.ticks.length=unit(-0.25, "cm"),
          #axis.title = element_blank(),
          #plot.subtitle = element_blank(),
          panel.border = element_rect(fill=NA, color="gray70", size=0.6),
          axis.ticks = element_line(color="gray70", size=0.3),
          plot.title = element_text(size=12, lineheight = 1, hjust=0, 
                                    margin = margin(3,0,0.5,0) ),
          plot.subtitle =  element_blank(),
          legend.position = "none",
          plot.margin = margin(0,0,1,1, "pt")
          
    )+
    labs(x="", y="", title=sub)+
    coord_cartesian( ylim=r, expand=0)#+  #
  #facet_wrap(~pc, nrow=1, scale="free_y")
  # if(left){
  #   g=g+theme(axis.text.y = element_text(size=9))
  # }
  # if(top){
  #   g=g+labs(subtitle=sub)+
  #     theme(plot.subtitle = element_text(size=13, lineheight = 1.2, hjust=0))
  # }
  # if(right){
  #   g=g+scale_y_continuous(sec.axis = dup_axis())+
  #     theme(
  #       axis.text.y = element_blank(),
  #       axis.text.y.right  = element_text(size=6, margin = margin(0,0,0,0)),
  #       axis.title.y.right = element_blank()
  #           )
  # }
  return(g)
}



do_median_timedist = function(xtest, nneighbors = 50){
  xtest = xtest %>% ungroup() %>% mutate(Date2=as.integer(Date2))
  difrange= (-(max(xtest$Date2)-min(xtest$Date2))): (max(xtest$Date2)-min(xtest$Date2))
  dyears=sort(unique(xtest$Date2))
  alldiffs = data.frame(dif=difrange, globaldif=0, neibdif_eucl=0, neibdif_cos=0, row.names =difrange) # gets incremented
  yeardiffs = matrix(0, nrow=length(difrange), ncol=length(dyears), dimnames =list(difrange, dyears))  # matrix of n diffs rows, y years, increment appropriate year
  # need to make sure that if using names for indexing, it's character not numeric
  pcm=xtest %>% select(starts_with("PC")) %>% as.matrix()
  xestim = data.frame(
    Date2 = xtest$Date2 ,
    expected_eucl = NA,
    actual_eucl = NA,
    estimate_eucl = NA,
    expected_cos = NA,
    actual_cos = NA,
    estimate_cos = NA
  )
  
  # record which works are same author, to make them na in calcs
  print(paste(Sys.time(), "doing sim matrices"))
  makena = vector("list", nrow(xtest))
  for(i in 1:nrow(xtest)){
    makena[[i]] = which(xtest$author_name==xtest$author_name[i])
  }
  # chunk sim&dist calculation, optimal for time&ram
  chunks = (1:nrow(xtest)) %>%  {split(., ceiling(seq_along(.)/5000))}
  neibs_eucl=vector("list", length(chunks))
  neibs_cos=neibs_eucl
  for(i in seq_along(chunks)){
    print(paste(Sys.time(), i))
    etmp = dist2(pcm, pcm[chunks[[i]],,drop=F], "euclidean","none")
    ctmp = sim2(pcm, pcm[chunks[[i]],,drop=F])
    for(j in seq_along(chunks[[i]]) ){ # make same author works na:
      etmp[makena[[ chunks[[i]][j] ]] , j] = NA
      ctmp[makena[[ chunks[[i]][j] ]] , j] = NA
    }
    neibs_eucl[[i]] = apply(etmp,2, function(x) head(order(x), nneighbors) ) %>% t
    neibs_cos[[i]] =  apply(ctmp,2, function(x) head(order(x, decreasing=T), nneighbors)) %>% t
    rm(ctmp,etmp); gc(verbose = F) # force memory management just in case
  }
  neibs_eucl = do.call(rbind, neibs_eucl)
  neibs_cos  = do.call(rbind, neibs_cos)
  
  
  print(paste(Sys.time(), "doing calcs"))
  for(i in 1:nrow(xtest)){
    if(i %% 10000 == 0){print(paste(Sys.time(), i))}
    # could be optimized to chunk the sim matrix calc and setting the NAs beforehand
    # i=which(xtest$target)[1] # debug
    xd = xtest$Date2[i]
    ad = xtest$Date2 - xd
    ad[makena[[i]]] = NA # if same author then na (from above)
    adt = table(ad)
    nd = names(adt)
    # update year global distribution; for each year, all distances (except same author)
    yeardiffs[nd, as.character(xd)] =  yeardiffs[nd, as.character(xd)] + adt # adds up to global total later
    # neibs_eucl = 
    #   (dist2(pcm[,,drop=F], pcm[i,,drop=F], "euclidean","none")[,1]) %>% 
    #   {.[is.na(ad)]=NA;.} %>% # disable self and same author
    #   {ad[order(.)]} %>%   # order time distances by compression space distance
    #   head(nneighbors)     # take top n
    etmp = ad[neibs_eucl[i,]]
    ctmp = ad[neibs_cos[i,]]
    
    # update actual of each work, 2 distance measures:
    xestim$actual_eucl[i] = median(etmp, na.rm=T)
    xestim$actual_cos[i] = median(ctmp, na.rm=T)
    # update counts in the distributions object:
    euct = table(etmp) 
    cost = table(ctmp)
    alldiffs[names(cost),  "neibdif_cos"] =  alldiffs[names(cost), "neibdif_cos"] + cost
    alldiffs[names(euct),  "neibdif_eucl"] =  alldiffs[names(euct), "neibdif_eucl"] + euct
  }
  gc(verbose = F)
  
  # calculate global expected probability distribution p for all dt's; each p(dt) = Dn_dt/D_dt
  alldiffs = alldiffs %>% 
    mutate(globaldif = rowSums(yeardiffs)) %>% 
    mutate(pdt_eucl = neibdif_eucl/globaldif) %>% 
    mutate(pdt_cos = neibdif_cos/globaldif)
  
  # get estimates of distances for each year
  yearestim_cos=setNames(rep(NA, length(dyears)), dyears)
  yearestim_euc=yearestim_cos
  difscores_euc = ((difrange * alldiffs$pdt_eucl)/ sum(alldiffs$pdt_eucl, na.rm=T) )*nneighbors 
  difscores_cos = ((difrange * alldiffs$pdt_cos)/ sum(alldiffs$pdt_cos, na.rm=T) )*nneighbors 
  for(i in seq_along(dyears)){
    if(any(yeardiffs[,i]>0)){
      # optimization: do weighted median, where year freqs are the weights; faster 
      yearestim_euc[i]  =  wtd.quantile( difscores_euc, yeardiffs[,i] , 0.5)   
      yearestim_cos[i]  =  wtd.quantile( difscores_cos, yeardiffs[,i] , 0.5)   
      # if expanding counts; same result
      #x = rep(difrange, yeardiffs[,i]); median(( (x*alldiffs[as.character(x), "pdt_cos"])/s)*nneighbors)
    }
    # the n_i_t: for each year, the counts; here can do lookup, given the year of target.
    # pdt and nit are vectors, need to make them match
    #Ex = ((pdt *  nit_year_y) / sum(pdt)) * n_neib
    #(table(allneibdiffs[[i]]$dif)*(as.numeric(names(table(allneibdiffs[[i]]$dif)))) ) %>% sum %>% {./100}
  }
  xestim =xestim %>%   # make sure character indices; otherwise will break silently
    mutate(expected_eucl = yearestim_euc[as.character(Date2)],
           expected_cos = yearestim_cos[as.character(Date2)] 
    ) %>% 
    mutate(
      estimate_cos =  actual_cos  - expected_cos,
      estimate_eucl = actual_eucl - expected_eucl
    )
  gc(verbose=F)
  return(xestim)
}


soldinsets = function(p,let, h=hdotssold){
  x=h %>% filter(pc==p) %>% 
    filter(value >=-20, value <= 10) %>% 
    mutate(dg = cut(Date2,90,include.lowest=T),pg= cut(value,60)) %>% 
    group_by(dg,pg) %>% 
    summarise(sold=mean(sold), count=n(), Date2=mean(Date2), value=mean(value)) %>% 
    ungroup() %>% 
    filter(count>1) %>% 
    mutate(count=case_when(count>=20 ~ 20L, T~count)) %>% 
    mutate(count=rescale(count,to=c(0.1,1)))
  
  ggplot(x, aes(dg ,pg ,fill=sold, alpha=count))+
    #geom_hline(yintercept = levels(x$pg)[seq(1,200, length.out=7)],color="gray70", size=0.1)+
    geom_tile(size=0)+
    annotate("text", x=-Inf, y=-Inf, hjust=-0.01,vjust=-0.2, label=paste0(let, " ", p, ", sold vs not sold"), size=4)+
    scale_fill_gradientn(limits=c(0,1),colors = diverging_hcl(palette="Blue-Red 3",n=11) %>% {.[6]="gray95";.})+
    scale_alpha(range=c(0.03,1))+
    coord_cartesian(expand=0) +
    theme_void()+
    theme(legend.position = "none", 
          panel.background = element_rect(color=NA,fill="white"),
          panel.border = element_rect(color="black", fill=NA))
}
