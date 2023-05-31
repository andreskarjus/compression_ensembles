############################### #
# Scripts to accompany the paper 
# "Compression ensembles quantify aesthetic complexity and the evolution of visual art"
# by Karjus et al 2023
# This software is provided as-is, without any implication of warranty.
# If you make use of this software or the precomputed datasets, kindly cite the paper.
#
# To start, download all the scripts and the precomputed objects of interests, and place
# them into a folder, and define the full path to this folder below.
# This script has 2 parts: using precomputed vectors as discussed in the paper,
# and a section for computing your own vectors.
############################## #



# Define parameters
params = list(
  # Define full path to the folder where you placed the scripts and the precomputed objects
  projectfolder = "C:/Users/Andres/korpused/artcorpora/art500kmodel"
  
  # Compression ensemble parameters - leave default to replicate paper:
  ,pixeltotal=160000
  ,pixelsize=400
  ,coefs = c(1, 0.4)
  ,extracoefs=c( 0.2, 0.1)
  ,minsize=0.5
  
  # The other folder paths are only needed if doing multicore processing
  , nfree = 1   # number of cores to not use/leave free if multicore   
  ,corpusfolder = ""  # folder with images, e.g. art500k
  ,modelfolder = ""   # folder where to save the computed ensemble vectors (as RData)
  ,chunksize=1000     # n images in each multicore run (lower it if mem issues)
)


# Load functions and load packages
# Will install packages if missing
source(file.path(params$projectfolder, "compression_ensembles_functions.R"))




#### Part 1: results in this paper, using precomputed embeddings ####


#### Precomputed art500k (wikiart) and hic et nunc vectors ####

load(file.path(params$projectfolder, "artmodel.RData")) # art500k incl their metadata
load(file.path(params$projectfolder,"hicmodel.RData"))  # hic et nunc
load(file.path(params$projectfolder, "artcol.RData"))   # precomputed avg color for each image



#### Figure 1 large UMAP ####

# Either calculate a new umap (takes a while; also, plot range probably need adjustment)
bigumap = umap(dat %>% select(starts_with("TR_")) %>% scale %>% as.data.frame())$layout %>% 
  cbind(dat %>% select(painting_name, Date2, Path2), . )
bigumap$artcol = artcol

# Or, load an already pre-computed umap:
load(file.path(params$projectfolder, "bigumap.RData"))
# The first two V1, V2 columns are the UMAP coordinates
# This is slightly zoomed in on the central cloud of data points.
ggplot(bigumap %>% filter(V2 > -5) , aes(V1, V2))+
  geom_point(size=0.2,alpha=1, color=bigumap$artcol[bigumap$V2>-5], shape=15, fill=NA)+
  scale_x_continuous(limits = c(-8,9), expand=c(0,0))+
  scale_y_continuous(limits = c(-5.5,6.4), expand=c(0,0))+
  theme_void()+
  theme(panel.background = element_blank(), 
        plot.background = element_rect(color="gray50", fill="gray50"),
        panel.border = element_blank(),
        plot.margin = margin(0.02,0.02,0.02,0.02))




#### Trends over time graph ####

# Run all of this to generate the results and graph

pcal = do_common_pca(dat, hdat)
attr(pcal, "topvars")[1:2]

dtart = do_mean_sd_space(pcal$art, miny=1500,maxy=2000,minwin=10,maxwin=25,target=1000, maxpc=3)
dthicall = do_mean_sd_space(pcal$hic,miny=1, maxy=175, minwin=10,maxwin=25,target=1000, maxpc=3)
dthicsold = do_mean_sd_space(pcal$hic %>% filter(!is.na(tradeprice)), miny=1, maxy=175, minwin=10,maxwin=25,target=1000, maxpc=3)

adots = pcal$art %>% mutate(col=artcol) %>% filter(Date2 %in% c(1500:2000)) %>% select(Date2,col,Genre, Style,Field, Path2, PC1:PC3) %>% pivot_longer(c(PC1:PC3), names_to="pc")

slabs = pcal$art %>% mutate(s=gsub("high |northern ", "",Style,ignore.case = T)) %>%  mutate(s=gsub("^([^,\\(/;]+)(,|\\(|/|;).*","\\1", s)) %>% group_by(s) %>% summarise(Date2=median(Date2), n=n()) %>% filter(Date2>1500) %>%  filter( (Date2>1900 & n>3000) | (Date2>1800 & n>1000) | (Date2<1800 & n>100) | (Date2 >1980 & n>300) ) %>% mutate(s=gsub(" ", "\n",s) %>% gsub("\n$","",.) %>% gsub("-","\n",.))  %>% filter(!is.na(s)) %>% arrange(Date2) %>% 
  filter(s!= "Post\nImpressionism") # not enough space
hdotsall = pcal$hic %>% mutate(col=hcol) %>% select(Date2,col,Path2, PC1:PC3) %>% pivot_longer(c(PC1:PC3), names_to="pc") 
hdotssold = pcal$hic %>% mutate(sold=case_when(is.na(tradeprice)~0,T~1)) %>% select(sold,Date2, PC1:PC3)  %>% pivot_longer(c(PC1:PC3), names_to="pc")

# do median image per day
x = pcal$hic %>% filter(!is.na(tradedate))
px = c()
for(i in sort(unique(x$Date2))){
  xx=x %>% filter(Date2==i)
  m= xx %>% select(starts_with("PC")) %>% apply(2,median) %>% matrix(nrow=1)
  px[i] = dist2(xx %>% select(starts_with("PC")) %>% as.matrix(),m, "euclidean", "none")[,1] %>% which.min() %>% {xx[.,"Path2"]}
}
pxl = lapply(c(1:25, 151:175), function(x){
  image_read(px[x]) %>% 
    image_resize(geometry_size_pixels(h=70))  %>% 
    image_annotate(x, size = 21, gravity = "northwest", color = "black") %>% 
    image_annotate(x, size = 20, gravity = "northwest", color = "white") %>% 
    image_border(color="white",geometry="3x3")
})
pxl2 = image_append(c(image_join(pxl[1:25]) %>% image_append(),
                      image_join(pxl[26:50]) %>% image_append() ), stack=T)

x = pcal$hic %>% mutate(col=hcol) %>%
  mutate(sold=case_when(is.na(tradeprice)~F,T~T)) %>% 
  filter(Date2 %in% (127:131), PC1 > -13, PC1 < -7) %>% 
  group_by(creatoralias) %>% mutate(creator=case_when(n()<10~"other",T~creatoralias))

g = 
  medplot(adots %>% filter(pc=="PC1"), dtart %>% filter(pc=="PC1"),r = c(-5,10), sub="(A) 68k artworks from art500k/Wikiart, 1500-2000. PC1 ~ texture & detail complexity") +
  geom_label(aes( x=1799, y=-1.1+0.7, label="Realism, Impressionism vs\nprevalence of Abstract art"), size=4.3, lineheight=1, vjust=1, hjust=1, color="black",inherit.aes = F, label.size = NA) +
  annotate("segment", x=1799,y=c(-1.3,-2)+0.5, xend=c(1875, 1975), yend=c(2, -0.1), alpha=1, size=0.3, color="black")+
  geom_text(aes(y=-4.89, label=s),data=slabs[seq(1,nrow(slabs),2),] ,angle=90,hjust=0,vjust=-0.1, color="black",lineheight=0.6, size=3.2, alpha=1)+
  geom_text(aes(y=9.89, label=s),data=slabs[seq(2,nrow(slabs),2),] ,angle=90,hjust=1, color="black", size=3.2,lineheight=0.6, alpha=1)+
  scale_x_continuous(breaks=seq(1500,1900,100)) +
  
  
  medplot(hdotsall %>% filter(pc=="PC1"), dthicall %>% filter(pc=="PC1"), r=c(-20,10), sub="(C) Hic et Nunc, March-August 2021, 50k objects. PC1 ~ texture & detail") +
  #
  geom_label(aes( x=40, y=8.9, label="complex digital art"), size=4.3, lineheight=1, vjust=1, hjust=0, color="black", label.size=NA)+
  annotate("segment", x=40,y=8.2, xend=c(32), yend=c(5), alpha=1, size=0.3, color="black")+
  #
  geom_label(aes( x=80, y=-19.8, label="CryptoPunks-like avatars\nseries e.g. AI Pokemon,\nDino Dudes, NFT-People"), size=4.3, lineheight=1, vjust=0, hjust=1, color="black", label.size=NA)+
  annotate("segment", x=81,y=-16, xend=c(109, 105), yend=c(-14.7,-8), alpha=1, size=0.3, color="black") +
  scale_y_continuous(breaks=seq(-20,10,5))+
  
  
  medplot(adots %>% filter(pc=="PC2"), dtart %>% filter(pc=="PC2"), r=c(-5,10), sub="(B) 68k artworks from art500k/Wikiart, 1500-2000. PC2 ~ overall compressibility" )+
  geom_label(aes( x=1750.5, y=9.5, label="Cubism, Expressionism, Surrealism, etc.\n vs Rococo portraits"), size=4.3, lineheight=1.2, vjust=1, hjust=1, color="black", label.size=NA)+
  annotate("segment", x=1750.6,y=c(8.1, 9.1), xend=c(1775, 1950), yend=c(-2, 7.5), alpha=1, size=0.3, color="black") +
  scale_x_continuous(breaks=seq(1500,1900,100))+
  
  medplot(hdotsall %>% filter(pc=="PC2"), dthicall %>% filter(pc=="PC2"), r=c(-20,10), sub="(D) Hic et Nunc, March-August 2021, 50k objects. PC2 ~ overall compressib.") + scale_y_continuous(breaks=seq(-20,10,5)) +
  inset_element(p=soldinsets("PC1", "(E)"), left = 0.01,right = 0.38,bottom = 0.01,top=0.35, align_to = "panel")+
  inset_element(p=soldinsets("PC2", "(F)"), left = 0.4,right = 0.779,bottom = 0.01,top=0.35, align_to = "panel") + plot_layout(ncol=2, widths=c(0.55,0.45)) -
 
  ggplot()+annotation_raster(pxl2,0,1,0,1)+labs(title="(G) Hic et Nunc median sold image, days 1-25 (top) and and 151-175 (bottom row)")+theme_void()+theme(plot.title = element_text(size=12,margin = margin(3 ,0,0,0)), plot.margin=margin(0,0,0,0))+coord_cartesian(expand=0)+
  
  plot_layout(nrow=2, heights = c(1,0.12*1.2)) + plot_annotation(theme=theme(plot.margin=margin(0,0,0,0)))


ggsave(file.path(params$projectfolder, "fig_dynamics.pdf"), g, height=8, width=13, device=cairo_pdf) 



#### Predict hic et nunc sales ####

datx = hdat %>%
  mutate(xvar=case_when( !is.na(tradeprice) ~ "sold",
                         T~"notsold"
  ) %>% as.factor() ) %>%
  mutate(Date2 = (hdat$timestamp %>% as.Date() %>% difftime(timestamp %>% as.Date() %>% min , units = "days") %>% as.integer) + 1) %>% 
  select(xvar, tradeprice, Date2, starts_with("TR_", ignore.case = F)) %>% group_by(xvar) %>% slice(sample(1:n(), pmin(n(), 25000)) ) %>% ungroup()
count(datx, xvar)
xhic = detection_accuracy(datx, xvars=colnames(datx)[-(1:3)], nclass=25000, ntest=5000,maxtrain=20000, nruns=500, trainseq = c(20000), kappa=F, debug=T, prefix = "o", doconfusion = F )

datx = datx %>% mutate(train= runif(n())<0.8  )
frm=paste( setdiff(colnames(datx),c("Date2", "train", "xvar")), "*Date2") %>% paste(collapse="+") %>% paste("xvar ~ ", .) %>% as.formula()
m = (glm(frm, data=datx %>% filter(train), family="binomial"))
(ifelse(predict(m, newdata = datx %>% filter(!train), "response")>0.5,yes= "sold",no= "notsold") == (datx %>% filter(!train) %>% pull(xvar))) %>% sum() %>% {./nrow(datx %>% filter(!train))}

pricedat = hdat %>% filter(!is.na(tradeprice)) %>% mutate(logprice=log(tradeprice))
pricedat = cbind(pricedat %>% select(!starts_with("TR_")), 
                 pricedat %>% select(starts_with("TR_")) %>% scale()
)
frm = grep("TR_", colnames(hdat), value=T) %>% paste(.,collapse=" + ") %>%  paste("logprice~", .)
frm= grep("TR_", colnames(hdat), value=T) %>% paste(.,collapse=" *Date2 + ") %>%  paste("logprice~", .) %>% as.formula()
summary(lm(frm, pricedat %>% filter(!is.infinite(logprice)) )) 



#### Temporal resemblance ####

xestim100 = do_median_timedist(dat, nneighbors = 100)
datx0 = cbind(dat, xestim100[,-1]) %>% filter(Date2 %in% 1800:1990) %>% 
  group_by(author_name) %>% mutate(Daten = rescale(Date2, 0:1))
datx0$res_eucl = gam(actual_eucl ~ s(Date2), data=datx0) %>% residuals()
datx0$res_cos = gam(actual_cos ~ s(Date2), data=datx0) %>% residuals()

# artist subset
artists = dat %>% 
  filter(!grepl("Papalucas|Hiroshige|Gustave Dore", author_name)) %>%  # patchy data
  group_by(author_name) %>% 
  filter(n()>=90, min(Date2)>=1800, 
         max(Date2)<=1990 , 
         abs(min(Date2)-max(Date2))<=60,
         abs(min(Date2)-max(Date2))>=30,
         n_distinct(Date2)>=12,  # leave out clustered (bad dating?) ones
  ) %>% #summarise(max(Date2-lag(Date2), na.rm=T))
  arrange(Date2) %>%  # for lag calc
  filter(max(Date2-lag(Date2), na.rm=T)<=10) %>% # to avoid gap ones 
  count %>% arrange(desc(n)); artists %>% nrow()
##examples ##
## calculate individual slopes
acoefs = list(); preds = tibble()
datx = list()
for(a in unique(artists$author_name)){
  d = datx0 %>% filter(author_name %in% a) %>% arrange(Daten)
  g = gam(res_cos~s(Daten, bs='cs',k=4, fx = T), data=d[-c(1, nrow(d)),]) 
  # cut first/last, often atypical, also gaps (misdating, slow in start and end?)
  d$resid = as.numeric(abs(d$res_cos - predict(g, d)))
  datx[[a]]=d
  acoefs[[a]] = coef(g)
  predseq=seq(0,1,length.out=1000)
  pr=predict(g, newdata=data.frame(
    Daten=predseq), se.fit=T)
  preds = rbind(preds, 
                tibble(
                  #s=seq(0,1,length.out=nrow(d))
                  p=pr$fit,
                  upr=pr$fit+(pr$se.fit*2),
                  lowr=pr$fit-(pr$se.fit*2),
                  se=pr$se.fit,
                  s=predseq, 
                  a=a, m=median(d$res_cos), sds=sd(d$res_cos), 
                  r2=summary(g)$r.sq, mind = min(d$Date2), maxd= max(d$Date2)
                )[30:969,]
  )
}; acoefs = do.call(rbind, acoefs); datx = do.call(rbind, datx)




#### Evals: human norms ####

# Load compression scripts and packages; installs pkgs if necessary
source(file.path(params$projectfolder, "artfcompressor_functions.R")) 

# Place the folders with the Multipic and Fractals images (see links in the paper)
# in the project dir, place their respective datasets in a directory that is also
# a subfolder of the project dir, and define their (relative) paths here (so just folder names)
evals = do_all_evals(projectfolder=params$projectfolder,
                     multipicdata = "",
                     multipicnorms = "" ,
                     fractaldata = "",
                     fractalnorms = ""
                     )



#### Evals: art classifier ####

# drawing vs painting
datx = dat %>% mutate(Field=tolower(Field)) %>%
  mutate(xvar=case_when( (grepl("^oil", Field) & !grepl(",",Field)) ~ "oil painting",
                         ( grepl("^drawing", Field) & !grepl(",",Field)) ~ "drawing",
                         T~NA_character_
  )) %>% dat
filter(Date2>=1600, !is.na(xvar)) %>% select(xvar, starts_with("TR"))  %>% filter(xvar %in% c("drawing") | (grepl("^oil",xvar) & !grepl(",",xvar)) ) %>% group_by(xvar) %>% slice(sample(1:n(), pmin(n(), 1500)) ) %>% ungroup()
count(datx, xvar)
vars_oil = approx_importance(datx, nclass=1100, ntest=100, nruns=1000)
xoil = detection_accuracy(datx, xvars=vars_oil, nclass=1100, ntest=100,maxtrain=1000, nruns=500,
                          trainseq = c(10,100,1000), kappa=F, debug=F, prefix = "o", doconfusion = F )

# style
datx = dat %>% group_by(author_name) %>% filter(n()>10) %>% ungroup() %>% rename(xvar=Style) %>% select(xvar, starts_with("TR")) %>% group_by(xvar) %>% filter(n()>=1100, !is.na(xvar)) %>% ungroup()
datx$xvar %>% unique() %>% length 
vars_style = approx_importance(datx, nclass=1100, ntest=100, nruns=1000)
xstyle = detection_accuracy(datx, xvars=vars_style, nclass=1100, ntest=100,maxtrain=1000, nruns=500, trainseq =c(10,100,1000), kappa=F, prefix = "s", debug=F, doconfusion=T )

# century; 
datx = dat %>% group_by(author_name) %>% filter(n()>10) %>% ungroup() %>%
  #mutate(xvar=cut(Date2, seq(1000,2050, 50), right = F, include.lowest = T, dig.lab = 4) )  %>%
  #mutate(xvar=as.factor(round(Date2/100)*100)) %>%
  mutate(xvar=as.factor(floor(Date2/100)*100)) %>%
  select(xvar, starts_with("TR")) %>% group_by(xvar) %>% filter(n()>=1100) %>%
  slice(sample(1:n(), pmin(n(), 3000)  )) %>% ungroup() %>% droplevels()
datx$xvar %>% unique() %>% length
nrow(datx)
vars_year = approx_importance(datx, nclass=1100, ntest=100, nruns=1000)
xyear = detection_accuracy(datx, xvars=vars_year, nclass=1100, ntest=100, maxtrain=1000, nruns=500,
                           trainseq =c(10,100,1000), kappa=F, debug=F, prefix = "y" )

# landscapes vs portraits
datx = dat %>% filter(Date2>=1700) %>% select(Genre, starts_with("TR")) %>% rename(xvar=Genre) %>% filter(xvar %in% c("landscape", "portrait")) %>% group_by(xvar) %>% slice(sample(1:n(), 5000) ) %>% ungroup()
count(datx, xvar)
vars_landscape = approx_importance(datx, nclass=1100, ntest=100, nruns=1000)
xland = detection_accuracy(datx, xvars=vars_landscape, nclass=1100, ntest=100,maxtrain=1000, nruns=1000,
                           trainseq = c(10,100,1000), kappa=F, debug=F, prefix = "l", doconfusion = F )

# author
datx = dat %>% select(author_name, starts_with("TR")) %>% rename(xvar=author_name) %>% 
  group_by(xvar) %>% filter(n()>=110 & n()<250 )
datx$xvar %>% unique() %>% length 
vars_author = approx_importance(datx, nclass=110, ntest=10, nruns=NULL, fullcomb=T)
xauthor = detection_accuracy(datx, xvars=vars_author, nclass=110, ntest=10, maxtrain=100, nruns=10000,
                             trainseq =c(10,100), kappa=F,  prefix = "a" )

xdats = rbind(
  xyear %>% mutate(xclass="(G) 7 centuries" )  ,
  xoil %>% mutate(xclass="(I) Drawing or oil painting" ),
  xland %>% mutate(xclass="(H) Landscape or portrait")  ,  
  xauthor %>% mutate(xclass="(E) 91 artists" ),
  xstyle %>% mutate(xclass="(F) 13 style periods" )
) %>% group_by(xclass) %>% mutate(classcol = acc/max(acc))
xdats$xclass = xdats$xclass %>% as.factor()
# Results:
xdats








#### Part 2: How to compute your own vectors ####

# Load compression scripts and packages
source( file.path(params$projectfolder, "artfcompressor_functions.R")) 


#### Simple use case: just one image ####

# Given a full path to an image file, this will calculate its compression ensemble vector, of the current implementation of an (opinionated) set of transformations. Put this into a loop to do multiple, or use the multicore batch processing option below.
vec = docompressions("full/path/to/image", params)



#### Multicore processing ####

# This likely needs some tweaking for particular datasets and hardware, if used
# Reqires an RData file meta_with_paths.RData in the project folder.
# The default path column in the dat file is Path2
# If running on a linux server, might need to run this to install packages from source:
# linuxinstaller()

if(dir.exists(params$projectfolder)){
  source(file.path(params$projectfolder, "artcompressor_functions.R"))
  load(file.path(params$projectfolder, "meta_with_paths.RData"))
  dat$Path2 = file.path(params$corpusfolder, dat$Pathb)
  chunks = split(1:nrow(dat), ceiling((1:nrow(dat))/params$chunksize) )
  comps=vector("list", length(chunks))
  for(i in seq_along(chunks) ){
    errorprinter(paste(Sys.time(), i ))
    tryCatch({
      comps[[i]]  = docompressions_multicore(dat[chunks[[i]],"Path2",drop=F], params)
      save(comps, file=file.path(params$modelfolder, "comps.RData"))
      gc(verbose = F, full = T)
    }, error=function(e)errorprinter(e), warning=function(e)errorprinter(e) )
  }
  
  #dat2 = do_pca(comps %>% unlist(recursive = F), dat)
  save(comps, file=file.path(params$projectfolder, "comps.RData") )
  errorprinter("done")
} else {
  stop("data folder missing")
}




#### Accessing art corpora ####

# This repository does not host any of the actual image files; we only provide the precomputed
# ensembles, and functions to calculate ensembles based on image datasets.
# The art500k subset (dubbed "Historical") can be accessed on their projet website.
# In principle, the metadata dataframe also includes the web links of the original images,
# so a corpus could possibly be populated by accessing those.
# The hic et nunc dataset was acquired via the now defunct version of the original marketplace.


##################################################### #
