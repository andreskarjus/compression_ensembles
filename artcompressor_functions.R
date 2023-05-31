############################### #
# Functions to accompany the script file for the paper 
# "Compression ensembles quantify aesthetic complexity and the evolution of visual art"
# by Karjus et al 2023
# This software is provided as-is, without any implication of warranty.
# If you make use of this software or the precomputed datasets, kindly cite the paper.
#
# This is the compressor functions file, intended to be sourced from the script file.
############################## #



#### load packages, installs if needed ####

p=c("magick", "dplyr", "entropy", "gifski", "ggplot2", "colorspace", "fractaldim", "reshape2", "parallel") 
install.packages(setdiff(p, rownames(installed.packages()))) 
print(sapply(p, require, character.only = TRUE))

linuxinstaller=function(){
  p=c("magick", "dplyr", "entropy", "gifski", "ggplot2", "colorspace", "fractaldim", "reshape2") 
  install.packages(p,quiet = T, type = "source")
  print(sapply(p, require, character.only = TRUE))
}


#### functions ####

count_colors <- function(x, n, forplot=T) {
  d <- image_data(x) %>%  # now quantization before call
    apply(2:3, paste, collapse = "") %>% 
    as.vector %>% table() %>%  as.data.frame()
  d[,1] <- paste("#", d[,1], sep="")
  if(forplot){
    d = d[order(d$Freq, decreasing = T),]
    d$Freq = sapply(d$Freq, function(x){round(x/ n[2]) * n[2]})
    d = matrix(rep(d[,1], d[,2])[1:(n[1]*n[2])], ncol=n[1], nrow=n[2] )
  }
  return(d)
}
scramble_colors=function(x){
  x = image_data(x) %>% apply(2:3, paste, collapse = "")
  set.seed(1) # all noisification replacements consistently the same
  x[] = sample(paste0("#", x))
  return(t(x))
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}



fdcalc = function(x, forplot=F){
  # calculates fractal dimensions with a couple of different windows
  # first transforms/quantizes to blacknwhite (technically just any 2 colors) to keep it fast 
  # (also colors would be treated indexical which is not great)
  # can also return plot if needed for illustration
  x2 = image_quantize(x, 2, dither = F, colorspace = "lab") %>% 
    image_raster() %>% 
    mutate(col=as.numeric(as.factor(col))-1) %>% 
    reshape2::dcast(x ~ y, value.var = "col") %>% as.matrix()
  xf1 = fd.estimate(data=x2, methods = "filter1", 
                    window.size = round(mean(dim(x2))/10),
                    step.size = round(mean(dim(x2))/40)  
  ) # rough, but much faster and gives about the same result as 20x2
  xf2 = fd.estimate(data=x2, methods = "filter1", 
                    window.size = round(mean(dim(x2))/3),
                    step.size = round(mean(dim(x2))/4) ) # no sliding, just ~rule-of-thirds cut
  x2l = image_canny(x, "0x1+10%+40%") %>% 
    image_raster() %>% 
    mutate(col=as.numeric(as.factor(col))-1) %>% 
    reshape2::dcast(x ~ y, value.var = "col") %>% as.matrix()
  
  xf3 = fd.estimate(data=x2l, methods = "filter1", 
                    window.size = round(mean(dim(x2))/5),
                    step.size = round(mean(dim(x2))/6) ) # no sliding, just ~rule-of-thirds cut 
  
  
  # make it so that larger numbers (now closer to 0) correspond to higher complexity
  fvalues = 
    c(fractaldim1 = 0-mean(xf1$fd[,,1], na.rm=T), 
      fractaldim2 = 0-mean(xf2$fd[,,1], na.rm=T),
      fractaldim3 = 0-mean(xf3$fd[,,1], na.rm=T)
  )
  
  if(forplot){
    fvalues =  # if illustrating, just report actual values for simplicity
      c(fractaldim1 = mean(xf1$fd[,,1], na.rm=T), 
        fractaldim2 = mean(xf2$fd[,,1], na.rm=T),
        fractaldim3 = mean(xf3$fd[,,1], na.rm=T)
      )
    
    xlist=list(xf1,xf2,xf3)
    fimlist=vector("list", 3)
    xinfo = image_info(x)
    mn = xinfo[2:3] %>% min
    for(i in seq_along(xlist)){
      xplot = image_graph(xinfo$width, xinfo$height, antialias = F, clip = F, 
                          res = as.numeric(gsub("[^0-9]*([0-9]+).*","\\1", xinfo$density)) )
      x3 = xlist[[i]]$fd[,,1] %>% reshape2::melt() %>% 
        ggplot(aes(Var1,Var2, fill=value))+
        geom_tile()+
        scale_fill_gradient(low="gray5", high="gray90", na.value = "pink", guide = "none")+
        scale_x_continuous(expand = c(0,0))+
        scale_y_continuous(expand = c(0,0), trans = "reverse")+
        theme_void()+
        #coord_fixed()+
        theme(panel.background=element_rect(fill = NA,colour = NA),
              panel.border=element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"),
              axis.title = element_blank()
        ) 
      print(x3)
      dev.off()  # saves image_graph 
      fimlist[[i]]=xplot; rm(xplot)
    }
    
    return(list(fimlist=fimlist, fvalues=fvalues))
  } else {
    return(fvalues)
  }
  
}


colorfulness = function(tmp){
  # input must be from image_raster(baseline)$col %>% hex2RGB() 
  # implements M3 measure from "Measuring Colorfulness in Natural Images"
  # they show the rgb one is basically as good if not better than the lab/chroma one
  # but grayscale images get 0 though regardless of number of shades of gray.
  tmp = tmp %>% coords()  # rgb
  rg = tmp[,1]-tmp[,2]    # r-g
  yb = 0.5*(tmp[,1]+tmp[,2]) - tmp[,3]
  return(
    sqrt(var(rg)+var(yb))+(0.3*sqrt(mean(rg)^2 + mean(yb)^2))
  )
}


colorfulness2 = function(hexlab){
  # implements M2 measure from Measuring Colourfulness in Natural images
  # the other one gets range nicely but fails on grayscale; this one works too but seems more about saturation, but also works on gray, so will just use both
  tmp = hexlab %>% coords()
  ch = sqrt(tmp[,2]^2 + tmp[,3]^2)
  return(
    sqrt(var(tmp[,2])+var(tmp[,3])) +(0.94*mean(ch))
  )
}



transformer_example = function(x, pipe=F){
  x = image_convert(x, format="rgb", depth = 8, antialias = F) %>% image_resize(geometry_size_pixels(width  = 200))
  
  xinfo = image_info(x)
  # imgray = image_convert(x, colorspace = "Gray")
  # imblur10 = image_blur(x, 10,10)
  # imdesp10 = image_despeckle(x,10)
  imdesp5 = image_despeckle(x,5) 
  imdesp2 = image_despeckle(x,2)
  #imcanny = imdesp2  %>% image_canny("0x1+10%+40%") 
  tmp = c(
    x
   # , image_resize(x, geometry_size_pixels(width  = 4, preserve_aspect = F)) %>% image_modulate(brightness = 10000)
    #, imblur10 = image_blur(x, 30)
    ,morph = image_scale(x, geometry_size_pixels(width = xinfo$width/10, height = xinfo$height/10)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F)) %>% image_modulate(brightness=150)
    ,imq = x %>% image_modulate(brightness=300) %>% image_quantize(2, dither = F, colorspace = "hcl")  %>% image_convert(colorspace = "Gray")
    #,lines_edge10_gray = imdesp2  %>% image_quantize( 10, dither = F, colorspace = "Gray") %>%   image_edge(10)
   , image_lat(x) %>% image_negate()
   #, image_flatten(c(x, image_channel(x, "luminance")), "divide") 
    ,lines_color_comp = 
      #image_compare(x, imdesp5, fuzz = 5) %>% image_fx( expression = "round(p)") %>% image_convert(colorspace = "Gray") %>% image_negate()
      imdesp2  %>% image_canny() %>% image_negate()
  ) 
  tmp = image_border(tmp, "gray99", "5x5")
  if(pipe){
    return(image_montage(tmp[2:5],geometry ="80x100+0+0", tile = "2x2", bg = "white"))
  }
  return(
    image_append(tmp)
  )
}


transformer = function(x, fillcol="green"){
  # some of the transformations make an image more complex (e.g. via noise, emboss), some make it simpler (e.g. washing out colors)
  xinfo = image_info(x)
  mn = xinfo[2:3] %>% min
  {
    x2 = image_graph(xinfo$width, xinfo$height, antialias = F, clip = F, 
                     res = as.numeric(gsub("[^0-9]*([0-9]+).*","\\1", xinfo$density)) )
    x3 = image_ggplot(x) + geom_polygon(aes(x,y),data=circleFun(c(xinfo$width/2,xinfo$height/2),(mn)-(mn*0.05),npoints = 50), fill=fillcol)
    print(x3)
    dev.off()
  }
  rm(x3)
  f=10 # fill coefficient
  # prepare transforms used for multiple operations:
  imgray = image_convert(x, colorspace = "Gray") # %>% image_convert(colorspace = "rgb")
  imblur10 = image_blur(x, 10,10)
  imdesp10 = image_despeckle(x,10)
  imdesp5 = image_despeckle(x,5) 
  imdesp2 = image_despeckle(x,2)
  imcanny = imdesp2  %>% image_canny("0x1+10%+40%")  # 
  
  return(
    list(
      #  x  # compress without transformation too  # now in sizer + different algos
      #, image_blur(x, 2, 2) #
      blur10 = imblur10 #
      , blur30 = image_blur(x, 30, 30) #
      #, blur10_lighten = image_composite(imblur10, x, "lighten") #
      #, image_noise(x, "Gaussian") #
      #, image_modulate(x, saturation = 25)  # halved saturation/semi-grayscale
      
      , colors_brightness = image_modulate(x, brightness = 300) # washes all lighter colors to white
      #, image_modulate(x, brightness = 10)
      #, image_equalize(x)
      , colors_acos = image_fx(x, expression = "acos(p)")
      , colors_sqrt = image_fx(x, expression = "sqrt(p)")
      #, image_fx(x, expression = "cos(p)") %>% image_negate()
      #, image_fx(x, expression = "atan(p)")
      #, image_fx(x, expression = "(p^2)")
      , colors_p10 = image_fx(x, expression = "(p^10)") # high cor with acos;  # these are all slow and similar so will use only some
      , colors_round = image_fx(x, expression = "round(p)") # rounding/simplifying of rgb & hcl channels:
      #, image_fx(x, expression = "round(p)", channel = "lightness")
      , colors_roundchroma = image_fx(x, expression = "round(p)", channel = "chroma")
      #, colors_roundhue = image_fx(x, expression = "round(p)", channel = "hue")  # very highly correlated with the chroma one so not needed really
      , colors_quantize_bw_dither = image_convert(x, type = 'Bilevel') # bw dithered
      , colors_quantize_bw = image_quantize(image_modulate(x, brightness = 150),2, dither = F, colorspace = "Gray")   # bw plain
      , colors_quantize3 = image_quantize(x,3, dither = F, colorspace = "rgb") # try rgb too
      , colors_quantize5 = image_quantize(x,5, dither = F, colorspace = "lab")
      #, colors_quantize10 = image_quantize(x,10, dither=F, colorspace = "lab") # correlates with 4
      , colors_saturate = image_modulate(x, saturation = 5000) # max all rgb vals
      , colors_grayscale = imgray 
      # , color_channel_hue = image_channel(x, "hue")     # hcl channel values as grayscales
      # , color_channel_chroma = image_channel(x, "chroma") # 0.99 cor with just grayscale
      , color_chroma_divide = image_flatten(c(x, image_channel(x, "chroma")), "divide") # image values divided by the chroma channel values
      , color_luminance_lighten = image_flatten(c(x, image_channel(imblur10, "luminance")), "linearlight") # looks burnt
      , color_luminance_divide = image_flatten(c(x, image_channel(x, "luminance")), "divide") # looks like xray
      , color_darken_intensity = image_flatten(c(x, image_negate(imdesp10)), "DarkenIntensity")
      #, imdesp5 %>%  image_convert(colorspace = "Gray") #
      #, imdesp2
      , colors_add2 = image_flatten(c(x, x), "add") %>%  image_despeckle(1) # add image values
      
      , morph_despecle10 = imdesp10    # watercolory effect; also removes noise
      , morph_pixelate10 = image_scale(x, geometry_size_pixels(width = xinfo$width/10, height = xinfo$height/10)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F))
      , morph_pixelate20 = image_scale(x, geometry_size_pixels(width = xinfo$width/20, height = xinfo$height/20)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F))
      , morph_add3_pixelate = image_flatten(c(x, x, x), "add")  %>% # image_despeckle(12)  %>% 
        image_scale(., geometry_size_pixels(width = xinfo$width/3, height = xinfo$height/3)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F)) # pretty groovy
      
      , morph_oilpaint = image_oilpaint(imdesp10, radius = 5) # color blobs basically
      #, image_oilpaint(x, radius = 15)
      , morph_squares = image_morphology(x, 'Dilate', "Square", iterations = 7)
      
      , lines_color_conv =  image_convolve(x, kernel = "7x7: 10 -5 -2 -1 -2 -5 -10 -5 0 3 4 3 0 -5 -2 3 6 7 6 3 -2 -1 4 7 8 7 4 -1 -2 3 6 7 6 3 -2 -5 0 3 4 3 0 -5 -10 -5 -2 -1 -2 -5 -10", iterations = 2)  # Sobel filter
      , lines_color_comp = image_compare(x, imdesp5, fuzz = 5) %>% image_fx( expression = "round(p)")  # another edgedet, via comparsion to blur
      #, image_edge(imdesp5, 1)
      , lines_edge1_color  = image_quantize(x, 5, dither = F, colorspace = "lab") %>% image_edge(1)
      , lines_edge2_color = image_edge(x,2)  
      , lines_edge_lat = image_lat(x) # edge detection via local adaptive thresholding
      , lines_cartoon = image_flatten(c(imdesp2, imcanny %>% image_negate()), "multiply") # cartoon
      , lines_edge5_gray =  image_edge(imgray,5) # only low correlation with lines_color_edge1
      , lines_edge10_gray = imdesp2  %>% image_quantize( 10, dither = F, colorspace = "Gray") %>%   image_edge(10)  # lines, but also kinda black ink effect
      #  , lines_bw_edge20 = image_edge(imgray,20)
      , lines_division_gray = image_flatten(c(imgray %>% image_blur(), imgray), "divide") # grayscale lines but also kinda paper cutout or chalk effect
      , lines_bw_canny = imcanny   # contours
      #, image_flatten(c(x, image_edge(imgray,10) ), "multiply")  # cartoony
      #, image_convolve(x,'DoG:0,0,2') %>% image_negate()  # similar to the edge 1+grayscale
      , lines_hough50 = image_hough_draw(imcanny,"50x50+70", color ="black", bg = "white") # composition angles, blank if no clear lines
      , lines_hough40 = image_hough_draw(imcanny,"40x40+20", color ="black", bg = "white")  # angles, more forgiving
      # "The Hough transform as it is universally used today was invented by Richard Duda and Peter Hart in 1972, who called it a "generalized Hough transform""
      #, image_quantize(x, 5, dither = F)  %>% image_hough_draw(color ="black", bg = "white") # angles
      #, image_convert(x, colorspace="gray") %>% image_threshold(type = "white", threshold = "50%") %>%image_threshold(type = "black", threshold = "50%") %>% image_hough_draw(color ="black", bg = "white")   # thought it'd be better but basically same
      # image_edge(image_convert(image_despeckle(x,5),colorspace = "Gray"),20)%>% image_negate() %>%  image_hough_draw(color ="black", bg = "white") # -> maybe best, for getting many angles, but slower 
      
      , flood_hole = x2   # central white hole
      , flood_thirds = image_fill(x, fillcol, point = geometry_point(xinfo$width/3, xinfo$height/3), f) %>% image_fill(fillcol, point = geometry_point( xinfo$width/1.5, xinfo$height/3),f) %>% image_fill(fillcol, point = geometry_point(xinfo$width/3, xinfo$height/1.5),f)  %>% image_fill(fillcol, point = geometry_point( xinfo$width/1.5, xinfo$height/1.5),f)  # fill rule of thirds crossings (has effect if something big on those points)
      , flood_corners = image_fill(x, fillcol, point = geometry_point(0,0), f) %>% image_fill(fillcol, point = geometry_point( xinfo$width-1, xinfo$height-1),f) %>% image_fill(fillcol, point = geometry_point( 0, xinfo$height-1),f)  %>% image_fill(fillcol, point = geometry_point( xinfo$width-1, 0),f)  # homogenize background/borders
      , flood_centre = image_fill(x, fillcol, point = geometry_point( xinfo$width/2, xinfo$height/2), f*2) # fill centre
      
      #, image_fuzzycmeans(x,smoothing = 1)  # different kinds of blobs
      #, image_fuzzycmeans(x,smoothing = 4) # too slow though
      , emboss_col1 = image_emboss(x, 1,0.1) # correlates with the other but not 100%. also just looks pretty.
      #, emboss_col4 = imdesp5 %>% image_emboss(4,1) # bumpmap looking thing; also kinda line/form detection (which it somewhat correlates with)
      , emboss_conv_grayd = image_convolve(imgray %>% image_negate() , 'DoG:0,0,2', 1,2, "20%")
      , emboss_conv_grayp = image_convolve(imgray, kernel = "Prewitt", 1, 1, "10%") 
      , emboss_gray4 = image_emboss(imgray, 4,1) # grayscale bumpmap thing
      , emboss_modulate = image_flatten(c(x, x), "modulate") %>% image_despeckle(1)
      #, emboss_multiply = image_flatten(c(x, image_quantize(x,3, dither = T, colorspace = "lab")), "multiply")
      
      , fx_deskew_zoom = image_deskew(x, threshold = 40) %>% image_crop(geometry = geometry_size_pixels(xinfo$width*0.7,  xinfo$height*0.7), gravity="center") %>% image_resize(geometry = geometry_size_pixels(xinfo$width,  xinfo$height)) # find position where lines most horizontal+vertical and recrop
      #, image_implode(x,-1) 
      , fx_implode = image_implode(x,1.01) %>% image_quantize(15, dither=F, colorspace = "lab") # weird central gravity well thing
      , fx_noise = image_noise(x, "Multiplicative")
      , fx_scramble = scramble_colors(x) %>%  image_read()  # scramble pixels (bc if already noise/abstract then size same)
      , fx_stripes = count_colors(image_quantize(x,20,dither=F), c(xinfo$width, xinfo$height)) %>% image_read() # get and sort (quantized) colors by frequency
    )
  )
}




transformer_small = function(x, fillcol=rgb(0,1,0)){
  # some of the transformations make an image more complex (e.g. via noise, emboss), some make it simpler (e.g. washing out colors), should this be recorded so can rotate the PCs accordingly to have the complex end always be positive, for ease of viewing? (no difference for other results)
  xinfo = image_info(x)
  mn = xinfo[2:3] %>% min
  {
    x2 = image_graph(xinfo$width, xinfo$height, antialias = F, clip = F,
                     res = as.numeric(gsub("[^0-9]*([0-9]+).*","\\1", xinfo$density)) )
    x3 = image_ggplot(x) + geom_polygon(aes(x,y),data=circleFun(c(xinfo$width/2,xinfo$height/2),(mn)-(mn*0.2),npoints = 50), fill=fillcol) # smaller circle in the smaller function
    print(x3)
    dev.off()
  }
  f=20 # fill coefficient          # larger fill here
  # prepare transforms used for multiple operations:
  imgray = image_convert(x, colorspace = "Gray")
  imblur10 = image_blur(x, 10,5)
  #imdesp10 = image_despeckle(x,10)
  #imdesp5 = image_despeckle(x,5) 
  imdesp2 = image_despeckle(x,2)
  imcanny = imdesp2  %>% image_canny("0x1+10%+40%")  # 
  fft =  image_fft(x) # fourier transform, 2 comps
  return(
    list(
      blur10 = imblur10 #
      # , blur30 = image_blur(x, 30, 40) #
      # , colors_brightness = image_modulate(x, brightness = 200) # washes all lighter colors to white
      , colors_acos = image_fx(x, expression = "acos(p)")
      # , colors_round = image_fx(x, expression = "round(p)") # rounding/simplifying of rgb & hcl channels:
      , colors_p10 = image_fx(x, expression = "(p^10)") 
      , colors_quantize_bw_dither = image_convert(x, type = 'Bilevel') # bw dithered
      , colors_quantize_bw = image_quantize(x,2, dither = F, colorspace = "Gray")   # bw plain
      , colors_quantize3 = image_quantize(x,3, dither = F, colorspace = "rgb") # try rgb too
      # , colors_quantize5 = image_quantize(x,5, dither = F, colorspace = "lab")
      , colors_grayscale = imgray 
      
      , morph_pixelate10 = image_scale(x, geometry_size_pixels(width = xinfo$width/10, height = xinfo$height/10)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F))
      , morph_pixelate20 = image_scale(x, geometry_size_pixels(width = xinfo$width/20, height = xinfo$height/20)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F)) 
      # , morph_oilpaint = image_oilpaint(imdesp2, radius = 2) # using 2 and 2 here
      
      , lines_edge1_color  = image_quantize(x, 10, dither = F, colorspace = "lab") %>% image_edge(1)
      , lines_cartoon = image_flatten(c(imdesp2, imcanny %>% image_negate()), "multiply") # cartoon
      , lines_edge10_gray = imdesp2  %>% image_quantize( 15, dither = F, colorspace = "Gray") %>%   image_edge(10)  # lines, but also kinda black ink effect
      , lines_edge5_gray = image_edge(imgray,5) 
      
      , lines_division_gray = image_flatten(c(imgray %>% image_blur(), imgray), "divide") # grayscale lines but also kinda paper cutout or chalk effect
      , lines_bw_canny = imcanny   # contours
      , lines_hough50 = image_hough_draw(imcanny,"50x50+70", color ="black", bg = "white") # composition angles, blank if no clear lines
      #, lines_hough40 = image_hough_draw(imcanny,"40x40+20", color ="black", bg = "white")  # angles, more forgiving
      
      , flood_hole = x2   # central white hole
      , flood_centre = image_fill(x, fillcol, point = geometry_point( xinfo$width/2, xinfo$height/2), f*3) # fill centre
      
      , emboss_gray4 = image_emboss(imgray, 4,-0.5) # grayscale bumpmap thing
      
      , fx_deskew_zoom = image_deskew(x, threshold = 50) %>% image_crop(geometry = geometry_size_pixels(xinfo$width*0.6,  xinfo$height*0.6), gravity="center") %>% image_resize(geometry = geometry_size_pixels(xinfo$width,  xinfo$height)) # find position where lines most horizontal+vertical and recrop
      # , fx_fourier1 = fft[1]
      # , fx_fourier2 = fft[2]
      # , fx_scramble = scramble_colors(x) %>%  image_read()  # scramble pixels (bc if already noise/abstract then size same)
    )
  )
}


transformer_bw = function(x, fillcol=rgb(0,1,0)){
  # assumes bw conversion has already been done
  xinfo = image_info(x)
  mn = xinfo[2:3] %>% min
  {
    x2 = image_graph(xinfo$width, xinfo$height, antialias = F, clip = F, 
                     res = as.numeric(gsub("[^0-9]*([0-9]+).*","\\1", xinfo$density)) )
    x3 = image_ggplot(x) + geom_polygon(aes(x,y),data=circleFun(c(xinfo$width/2,xinfo$height/2),(mn)-(mn*0.05),npoints = 50), fill=fillcol)
    print(x3)
    dev.off()
  }
  f=10 # fill coefficient
  # prepare transforms used for multiple operations:
  imblur10 = image_blur(x, 10,10)
  imdesp10 = image_despeckle(x,10)
  imdesp5 = image_despeckle(x,5) 
  imdesp2 = image_despeckle(x,2)
  imcanny = imdesp2  %>% image_canny("0x1+10%+40%")  # 
  fft =  image_fft(x) # fourier transform, 2 comps
  return(
    list(
      #  x  # compress without transformation too  # now in sizer + different algos
      #, image_blur(x, 2, 2) #
      blur5 =  image_blur(x, 5, 2)
      , blur10 = imblur10 #
      , blur30 = image_blur(x, 30, 30) #
      , desp5 = imdesp5
      #, blur10_lighten = image_composite(imblur10, x, "lighten") #
      #, image_noise(x, "Gaussian") #
      #, image_modulate(x, saturation = 25)  # halved saturation/semi-grayscale
      
      , colors_brightness = image_modulate(x, brightness = 300) # washes all lighter colors to white
      #, image_modulate(x, brightness = 10)
      #, image_equalize(x)
      #, colors_acos = image_fx(x, expression = "acos(p)")
      #, colors_sqrt = image_fx(x, expression = "sqrt(p)")
      #, image_fx(x, expression = "cos(p)") %>% image_negate()
      #, image_fx(x, expression = "atan(p)")
      #, image_fx(x, expression = "(p^2)")
      #, colors_p10 = image_fx(x, expression = "(p^10)") # high cor with acos;  # these are all slow and similar so will use only some
      #, colors_round = image_fx(x, expression = "round(p)") # rounding/simplifying of rgb & hcl channels:
      #, image_fx(x, expression = "round(p)", channel = "lightness")
      #, colors_roundchroma = image_fx(x, expression = "round(p)", channel = "chroma")
      #, colors_roundhue = image_fx(x, expression = "round(p)", channel = "hue")  # very highly correlated with the chroma one so not needed really
      # , colors_quantize_bw_dither = image_convert(x, type = 'Bilevel') # bw dithered
      #, colors_quantize_bw = image_quantize(image_modulate(x, brightness = 150),2, dither = F, colorspace = "Gray")   # bw plain
      #, colors_quantize3 = image_quantize(x,3, dither = F, colorspace = "rgb") # try rgb too
      #, colors_quantize5 = image_quantize(x,5, dither = F, colorspace = "lab")
      #, colors_quantize10 = image_quantize(x,10, dither=F, colorspace = "lab") # correlates with 4
      #, colors_saturate = image_modulate(x, saturation = 5000) # max all rgb vals
      #, colors_grayscale = imgray 
      # , color_channel_hue = image_channel(x, "hue")     # hcl channel values as grayscales
      # , color_channel_chroma = image_channel(x, "chroma") # 0.99 cor with just grayscale
      # , color_chroma_divide = image_flatten(c(x, image_channel(x, "chroma")), "divide") # image values divided by the chroma channel values
      , color_luminance_lighten = image_flatten(c(x, image_channel(imblur10, "luminance")), "linearlight") # looks burnt
      # , color_luminance_divide = image_flatten(c(x, image_channel(x, "luminance")), "divide") # looks like xray
      , color_darken_intensity = image_flatten(c(x, image_negate(imdesp10)), "DarkenIntensity")
      #, imdesp5 %>%  image_convert(colorspace = "Gray") #
      #, imdesp2
      #, colors_add2 = image_flatten(c(x, x), "add") %>%  image_despeckle(1) # add image values
      
      , morph_despecle10 = imdesp10    # watercolory effect; also removes noise
      , morph_pixelate10 = image_scale(x, geometry_size_pixels(width = xinfo$width/10, height = xinfo$height/10)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F))
      , morph_pixelate20 = image_scale(x, geometry_size_pixels(width = xinfo$width/20, height = xinfo$height/20)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F))
      #, morph_add3_pixelate = image_flatten(c(x, x, x), "add")  %>% # image_despeckle(12)  %>% image_scale(., geometry_size_pixels(width = xinfo$width/3, height = xinfo$height/3)) %>% image_scale(geometry_size_pixels(width = xinfo$width, height = xinfo$height, preserve_aspect = F)) # pretty groovy
      
      , morph_oilpaint = image_oilpaint(imdesp10, radius = 5) # color blobs basically
      #, image_oilpaint(x, radius = 15)
      , morph_squares = image_morphology(x, 'Dilate', "Square", iterations = 7)
      
      , lines_color_conv =  image_convolve(x, kernel = "7x7: 10 -5 -2 -1 -2 -5 -10 -5 0 3 4 3 0 -5 -2 3 6 7 6 3 -2 -1 4 7 8 7 4 -1 -2 3 6 7 6 3 -2 -5 0 3 4 3 0 -5 -10 -5 -2 -1 -2 -5 -10", iterations = 2)  # Sobel filter
      , lines_color_comp = image_compare(x, imdesp5, fuzz = 5) %>% image_fx( expression = "round(p)")  # another edgedet, via comparsion to blur
      #, image_edge(imdesp5, 1)
      #, lines_edge1_color  = image_quantize(x, 5, dither = F, colorspace = "lab") %>% image_edge(1)
      , lines_edge2_color = image_edge(x,2)  
      , lines_edge_lat = image_lat(x) # edge detection via local adaptive thresholding
      #, lines_cartoon = image_flatten(c(imdesp2, imcanny %>% image_negate()), "multiply") # cartoon
      , lines_edge5_gray =  image_edge(x,5) # only low correlation with lines_color_edge1
      , lines_edge10_gray = imdesp2  %>% image_quantize( 10, dither = F, colorspace = "Gray") %>%   image_edge(10)  # lines, but also kinda black ink effect
      #  , lines_bw_edge20 = image_edge(imgray,20)
      #, lines_division_gray = image_flatten(c(imgray %>% image_blur(), imgray), "divide") # grayscale lines but also kinda paper cutout or chalk effect
      , lines_bw_canny = imcanny   # contours
      #, image_flatten(c(x, image_edge(imgray,10) ), "multiply")  # cartoony
      #, image_convolve(x,'DoG:0,0,2') %>% image_negate()  # similar to the edge 1+grayscale
      , lines_hough50 = image_hough_draw(imcanny,"50x50+70", color ="black", bg = "white") # composition angles, blank if no clear lines
      , lines_hough40 = image_hough_draw(imcanny,"40x40+20", color ="black", bg = "white")  # angles, more forgiving
      # "The Hough transform as it is universally used today was invented by Richard Duda and Peter Hart in 1972, who called it a "generalized Hough transform""
      #, image_quantize(x, 5, dither = F)  %>% image_hough_draw(color ="black", bg = "white") # angles
      #, image_convert(x, colorspace="gray") %>% image_threshold(type = "white", threshold = "50%") %>%image_threshold(type = "black", threshold = "50%") %>% image_hough_draw(color ="black", bg = "white")   # thought it'd be better but basically same
      # image_edge(image_convert(image_despeckle(x,5),colorspace = "Gray"),20)%>% image_negate() %>%  image_hough_draw(color ="black", bg = "white") # -> maybe best, for getting many angles, but slower 
      
      , flood_hole = x2   # central white hole
      , flood_thirds = image_fill(x, fillcol, point = geometry_point(xinfo$width/3, xinfo$height/3), f) %>% image_fill(fillcol, point = geometry_point( xinfo$width/1.5, xinfo$height/3),f) %>% image_fill(fillcol, point = geometry_point(xinfo$width/3, xinfo$height/1.5),f)  %>% image_fill(fillcol, point = geometry_point( xinfo$width/1.5, xinfo$height/1.5),f)  # fill rule of thirds crossings (has effect if something big on those points)
      , flood_corners = image_fill(x, fillcol, point = geometry_point(0,0), f) %>% image_fill(fillcol, point = geometry_point( xinfo$width-1, xinfo$height-1),f) %>% image_fill(fillcol, point = geometry_point( 0, xinfo$height-1),f)  %>% image_fill(fillcol, point = geometry_point( xinfo$width-1, 0),f)  # homogenize background/borders
      , flood_centre = image_fill(x, fillcol, point = geometry_point( xinfo$width/2, xinfo$height/2), f*2) # fill centre
      
      #, image_fuzzycmeans(x,smoothing = 1)  # different kinds of blobs
      #, image_fuzzycmeans(x,smoothing = 4) # too slow though
      #, emboss_col1 = image_emboss(x, 1,0.1) # correlates with the other but not 100%. also just looks pretty.
      #, emboss_col4 = imdesp5 %>% image_emboss(4,1) # bumpmap looking thing; also kinda line/form detection (which it somewhat correlates with)
      #, emboss_conv_grayd = image_convolve(imgray %>% image_negate() , 'DoG:0,0,2', 1,2, "20%")
      , emboss_conv_grayp = image_convolve(x, kernel = "Prewitt", 1, 1, "10%") 
      , emboss_gray4 = image_emboss(x, 4,1) # grayscale bumpmap thing
      #, emboss_modulate = image_flatten(c(x, x), "modulate") %>% image_despeckle(1)
      #, emboss_multiply = image_flatten(c(x, image_quantize(x,3, dither = T, colorspace = "lab")), "multiply")
      
      , fx_deskew_zoom = image_deskew(x, threshold = 40) %>% image_crop(geometry = geometry_size_pixels(xinfo$width*0.7,  xinfo$height*0.7), gravity="center") %>% image_resize(geometry = geometry_size_pixels(xinfo$width,  xinfo$height)) # find position where lines most horizontal+vertical and recrop
      #, image_implode(x,-1) 
      , fx_implode = image_implode(x,1.01) %>% image_quantize(15, dither=F, colorspace = "lab") # weird central gravity well thing
      , fx_noise = image_noise(x, "Multiplicative")
      , fx_fourier1 = fft[1]
      , fx_fourier2 = fft[2]
      #, fx_scramble = scramble_colors(x) %>%  image_read()  # scramble pixels (bc if already noise/abstract then size same)
      #, fx_stripes = count_colors(image_quantize(x,20,dither=F), c(xinfo$width, xinfo$height)) %>% image_read() # get and sort (quantized) colors by frequency
    )
  )
}




sizer = function(y, algo="gif", quality=NULL, depth=8, compression=NULL){
  # gif ignores depth arg, only tiff uses the compression arg, only jpeg uses quality arg
  s = colMeans(
    rbind(
      sapply(y, function(x) image_write(x, raw(), format = algo, 
                                        quality = quality, depth=depth,compression=compression) %>% length ),
      sapply(y, function(x) image_write(image_rotate(x, 90), raw(), format = algo, 
                                        quality = quality, depth=depth,compression=compression) %>% length )
    )
  )
  names(s) = names(y)
  return(s)
  # png and gif correlate fairly well, but gif makes more sense: 
  # a noisified image of the same pixels takes more space than an organized one.
  # with gif, need to rotate and take mean, as some variation in those (in png&jpeg too but less)
}

fourier_sizer = function(x){
  # fft uses image width for both w and h of transformed square
  # so rotate here, send to sizer, take average
  fft  = image_fft(x) # fourier transform, 2 comps
  fftr = image_rotate(x, 90) %>% image_fft()
  return(
    c(
      (sizer(fft[1], "gif") +  sizer(fftr[1], "gif"))/2,
      (sizer(fft[2], "gif") +  sizer(fftr[2], "gif"))/2
    ))
}
file.path("C:/Users/Andres/korpused/artcorpora/art500k/")
# x = image_read("C:/Users/Andres/korpused/artcorpora/art500k/Artists2_Albert_Bierstadt_Autumn_Woods_1886.jpg") %>% image_resize(geometry_size_pixels(200))  # debugger
# x = image_read("C:/Users/Andres/korpused/artcorpora/art500k/Artists1_Abbott_Handerson_Thayer_Winter_Landscape__yAF_Q8QGlViTDw.jpg") %>% image_resize(geometry_size_pixels(200))  # debugger


# microbenchmark(image_read(c(fp,fp,fp,fp,fp,fp,fp,fp,fp,fp)), sapply(1:10, function(x)image_read(fp)), times = 10)
# microbenchmark(image_convert(c(x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x),colorspace = "Gray"), sapply(1:20, function(y)image_convert(x,colorspace = "Gray")), times = 10) # basically same, would buy couple of seconds on 60k
# microbenchmark::microbenchmark(image_write(x, raw(), "gif" ), image_write(x, raw(), "png" ), times = 100)  # png is faster but would win like 2h tops on the full 60k


exampleplot_docompressions = function(fp,params, ifolder="C:/Users/Andres/korpused/artcorpora/illustrate"){
  # fp="C:/Users/Andres/korpused/artcorpora/chase.jpg"
  pixeltotal=params$pixeltotal
  pixelsize =params$pixelsize
  coefs = params$coefs
  extracoefs = params$extracoefs
  minsize = params$minsize
  
  im = image_read(fp) %>% image_convert(format="rgb", depth = 8, antialias = F)
  xy = as.numeric(image_info(im)[2:3])
  
  # calculate best x*y ratio to get equivalent total px for baseline and then downsample
  newsize = list()
  for(i in seq_along(coefs)){
    if(coefs[i]==1){
      newsize[[i]]=round(xy * sqrt(pixeltotal / (xy[1] * xy[2])))
    } else {
      newsize[[i]] =  round(xy * sqrt( (pixeltotal*(coefs[i]^2)) / (xy[1] * xy[2])))
    }
  }
  newsize = newsize %>% lapply(., function(x) if(all(x <= xy)){x}else{xy}   )  # but make sure no upsampling
  
  if(all(newsize[[1]] == xy ) ){
    # if already exactly correct size, then skip the unnecessary resize step
    baseline = im %>%  image_flatten() 
  } else {
    baseline = image_resize(im, paste0(newsize[[1]][1],"x", newsize[[1]][2], "!") ) %>% 
      image_flatten() # flattening by itself to remove transparency
  }
  image_destroy(im)
  
  # get col triplets; will reuse multiple places below
  imcol = image_raster(baseline)$col
  
  fillcol = "white" # for example
  #   c(setdiff(
  #   tolower(rgb(c(0,1,1,0,0,1,0,1),c(1,0,1,0,0,0,1,1),c(0,1,1,0,1,0,1,0))),
  #   table(imcol) %>% {./prod(newsize[[i]])} %>% .[.>0.001] %>% names()
  # ), "white")[1]
  # 
  
  
  # for perceptual lightness/contrast:
  imhex = imcol %>% hex2RGB() 
  rm(imcol) 
  hexlab = imhex %>% as("LAB")
  lb =  hexlab %>%  coords() %>% .[,1]    # extract L channel
  
  # color range calc (see function for details)
  colorrange = colorfulness(imhex)
  rm(imhex)
  colorrange2 = colorfulness2(hexlab)
  rm(hexlab)
  
  contrastsd =  sd(lb/100) 
  contrastrange = table(round(lb)) %>%
    .[.>params$pixeltotal*0.001] %>% # exclude outliers (few super bright or dark pixels)
    names() %>% as.numeric() %>% range() %>% 
    {abs(.[2]-.[1])/100} # standardize to 0-1
  rm(lb)

  sizes = vector("list", 1) #; names(sizes)=as.character(coefs)
  bmpsize  = image_write(baseline, raw(), "rgb", depth=16) %>% length()
  baseinisize = bmpsize # will do other compressions in the end which will use this
  tmp = transformer(baseline, fillcol=fillcol)

    
  # smaller = more compressible, simpler;
  gftmp=sizer(baseline, "gif")
  sgif = sizer(tmp, "gif") / gftmp
  #spng = sizer(tmp, "png", depth=16) / pngtmp
    
  # compressions without transforms but different compressors
  notr = c(
    compress_gif=gftmp
  ) / baseinisize   # normalize by bitmap/uncompressed size
  
  sizes = c(
    notr,
    sgif 
  )
    
    # do fourier transforms - needs to be resized separately; will do only for baseline
      # only do with baseline size
      forms = c("colors_quantize3", "blur10", "lines_division_gray")
      fftsizes = c(
        fourier_sizer(baseline),
        fourier_sizer(tmp[[forms[1]]]),
        fourier_sizer(tmp[[forms[2]]]),
        fourier_sizer(tmp[[forms[3]]])
      ) / baseinisize
      names(fftsizes) = 
        c("fft1", "fft2",paste0(c("fft1", "fft2"), "_", rep(forms,each=2)))
    
      fftimgs=list()
      fftimgs[[1]] = image_fft(baseline)[1] %>% image_modulate(brightness=50000)
      fftimgs[[2]] = image_fft(baseline)[2] 
      ii=3
      for(i in 1:3){
        xf=image_fft(tmp[[forms[i]]])
        fftimgs[[ii]] = xf[1] %>% image_modulate(brightness=50000)
        ii=ii+1
        fftimgs[[ii]] = xf[2]
        ii=ii+1
      }
      names(fftimgs) = 
        c("fft1", "fft2",paste0(c("fft1", "fft2"), "_", rep(forms,each=2)))
  
  sizes2 = unlist(sizes, F, use.names = T) # to vector
  # standardized unique names now
  

  # add summary color distribution stats 
  colfreq = count_colors(image_quantize(baseline,200, dither=F, colorspace = "lab"), pixelsize, forplot=F) %>% {.[,2]/sum(.[,2])}

  #add entropy of estimated composition lines
  a=image_despeckle(baseline,2)  %>% image_canny()
  anums = a %>%
    image_hough_txt("40x40+30")  %>% strsplit("#") %>% unlist(F,F) %>%
    .[-(1:4)] %>% gsub("^ ([0-9]+) .*", "\\1",.) %>% as.numeric() %>%
    .[!is.na(.)]
  angleentropy = anums %>% entropy()
  # if no lines, returns 0 entropy (noise should yield noisy canny lines > many hough lines > high entropy; but not impossible this fails sometimes, if the canny algo  doesn't pick up any structure at all in the noise)

  angleplot = a %>%
    image_hough_draw("40x40+30", color = "black") %>%
    image_annotate("S(...)", gravity = "northwest", location="+20+10", size=30,kerning = 1.1, color="black", strokecolor = "black")
  
  
  fd = fdcalc(baseline, forplot = T)  # list of 2: list of 3 imgs; 1 vector
  
  # imblank = image_blank(xy[1], xy[2], color="white")
  # formulaplot = function(ifile, b=imblank){
  #   image_composite(b,
  #                   image_read(file.path(ifolder, ifile)) %>% 
  #                     image_resize(geometry=geometry_size_pixels(w=xy[1])))
  # }
  
  
  
  # adjust some for visibility - for example graphs only
  tmp$lines_bw_canny = tmp$lines_bw_canny %>% image_negate()
  tmp$lines_division_gray = tmp$lines_division_gray %>% image_contrast(0.9)
  # these both preserve absolute white and black but modulate midtones:
  tmp$colors_p10 = tmp$colors_p10 %>% image_modulate(brightness = 3000)
  tmp$colors_acos = tmp$colors_acos %>% image_negate() %>% image_modulate(brightness = 300) %>% image_negate()
  
  imlist = c(
    baseline
    , tmp
    , fftimgs
    , fd$fimlist[1:2] # exclude 3rd
    # , angleplot
    # , formulaplot("m2.png")
    # , formulaplot("m3.png")
  ) %>% lapply(function(x) x %>% 
                 image_convert(colorspace="RGB") %>% 
                 image_resize("x200") %>% 
                 image_border( "gray70", "1x1")# %>% image_border("white", "7x0") 
               )
  
  
  sizevec = c(
    sizes2        # transforms;  has names
    , fftsizes    # fourier on baseline & 3 transforms
    , fd$fvalues[1:2] # exclude 3rd, boring
  )
  if(length(sizevec) != length(imlist)){stop("mismatch")}
  sdcolfreq=sd(colfreq) 
  sdcolfreq=ifelse(is.na(sdcolfreq), 0, sdcolfreq)
  statvec = c(
      stats_angleentropy = angleentropy
    , stats_colorfulness_rgb = colorrange
    , stats_colorfulness_lab = colorrange2
    , stats_contrastsd = contrastsd
    , stats_contrastrange = contrastrange

    # all these inverse (1-) in implementation, but easier to put raw on plot
    , stats_colfreq_mean = mean(colfreq)
    , stats_colfreq_median = median(colfreq) 
    , stats_colfreq_sd = sdcolfreq
    , stats_colfreq_max = max(colfreq)
    , stats_colfreq_entropy = entropy::entropy(colfreq)
  )
  stattext = 
    paste0(c(
      "S(angles)    = ",
      "color_m1     = ",
      "color_m2     = ",
      "sd(contrast) = ",
      "range(contr) = ",
      "mean(cols)   = ",
      "median(cols) = ",
      "sd(cols)     = ",
      "max(cols)    = ",
      "S(cols)      = "
          ), round(statvec,2)) %>% paste0(collapse="\n")
  
  
  abr = c(", b=", rep(", c=", length(sizevec)-3), rep(", f=",2))
  glist=vector("list", length(imlist))
  for(i in seq_along(imlist)){
    glist[[i]] =
      ggplot()+
      annotation_raster(as.raster(imlist[[i]]),0,1,0,1)+
      theme_void()+coord_cartesian(expand=F)+
      labs(title=paste0(" ", names(sizevec[i]) %>% 
                          {ifelse(nchar(.)>11, paste0(substr(.,1,11),"."),.) },
                        abr[i], round(sizevec[i],2)))+ 
      theme(plot.title = element_text(size=5, margin=margin(1,0,0, 1)),
            plot.margin = margin(1,3,1,3))
  }
  
  glist = c(glist, list(
    # stats here
    ggplot()+
      annotate("text",size=1.7, hjust=-0.01, vjust=1.03, x=-Inf, y=Inf, label=stattext, lineheight=0.8, family="mono")+
      labs(title="statistics")+ 
      theme_void()+coord_cartesian(expand=F, xlim = c(0,1), ylim=c(0,1))+
      theme(plot.title = element_text(size=5, margin=margin(1,0,0, 1)),
            plot.margin = margin(0,0,0,0),
            panel.border = element_rect(fill="transparent", color="gray70", size=0.1))
    
    
  ))
  
  
  
  gg = wrap_plots(glist, ncol = 9)+plot_annotation(theme = theme(plot.margin = margin(0,0,0,0,unit="pt")))
  # ggsave("fig_illustrate.pdf", gg, width=9, height=7)
  attr(gg, "clabs") = c(
      names(tmp)
    , names(fftimgs)
    , names(fd$fimlist[1:2])
  )
  attr(gg, "namelabs")=c(names(sizevec),"stats" )
  return(gg)
  
}









# called by docompressions_multicore
docompressions = function(fp, params){
  print(paste(Sys.time(), fp))
  compvec=NA # if pipeline fails at any point, will return 1-length NA vector, easy to catch later
  tryCatch({
    # pixeltotal=250000; pixelsize=500
    pixeltotal=params$pixeltotal
    pixelsize =params$pixelsize
    coefs = params$coefs
    extracoefs = params$extracoefs
    #normsize = params$normsize
    #mainalgo = params$mainalgo
    minsize = params$minsize
    # if(is.null(normsize)){normsize="gif"} # if not specified, assume normalizer is baseline size
    # if(is.null(minsize)){minsize=0.5} # if not specified, assume size threshold is half of baseline
    # if(is.null(mainalgo)){mainalgo="gif"}  # by default, gif
    
    converttogray = params$gray
    if(is.null(converttogray)){converttogray=F}
    
    # read and force rgb 256 colors (like gif) to enforce same/similar behaviour regardless of backend and os
    im = image_read(fp) %>% 
      image_convert(format="rgb", depth = 8, antialias = F) %>% 
      image_quantize(256, dither = F, colorspace = "rgb")
    # depth arg doesn't actually do anything?
    
    xy = as.numeric(image_info(im)[2:3])
    
    # if input image too small, fail with a message, caught by trycatch; vector will be NA
    if(prod(xy) < (pixeltotal*minsize) ){
      stop(paste0("Image too small:",fp, "  ", xy[1],"x", xy[2]) )
    }
    
    # if image is intended to be grayscale, but might not be
    if(converttogray){
      im = image_convert(im, colorspace = "Gray")
    }
    
    # calculate best x*y ratio to get equivalent total px for baseline and then downsample
    newsize = list()
    for(i in seq_along(coefs)){
      if(coefs[i]==1){
        # newsize[[i]]=round(xy * sqrt(pixeltotal / (xy[i] * xy[i]))) # old, bug
        newsize[[i]]=round(xy * sqrt(pixeltotal / (xy[1] * xy[2])))
      } else {
        newsize[[i]] =  round(xy * sqrt( (pixeltotal*(coefs[i]^2)) / (xy[1] * xy[2])))
      }
    }
    newsize = newsize %>% lapply(., function(x) if(all(x <= xy)){x}else{xy}   )  # but makes sure no upsampling
    
    # print(paste(Sys.time(), "new sizes calculated"))
    
    if(length(newsize)>1){
      # second safeguard: will fail if file too small to be downsampled comparably
      if(all(newsize[[1]] <= newsize[[2]])){
        warning(paste(fp, "too small native resolution"))
        return(NA) # returns single NA
      }
    }
    # should check lengths and set row to be NAs after the parLapply part
    
    if(all(newsize[[1]] == xy ) ){
      # if already exactly correct size, then skip the unnecessary resize step
      baseline = im %>%  image_flatten() 
    } else {
      baseline = image_resize(im, paste0(newsize[[1]][1],"x", newsize[[1]][2], "!") ) %>% 
        image_flatten() # flattening by itself to remove transparency
    }
    image_destroy(im)
    
    # get col triplets; will reuse multiple places below
    imcol = image_raster(baseline)$col
    
    # infer suitable color for flood-filling, so that it doesn't overlap with the image's primary colors; neon green or magenta probably least likely; but include backup white if none fit (might be some very neon/primary color image which would have all these)
    fillcol = c(setdiff(
      tolower(rgb(c(0,1,1,0,0,1,0,1),c(1,0,1,0,0,0,1,1),c(0,1,1,0,1,0,1,0))),
      table(imcol) %>% {./prod(newsize[[i]])} %>% .[.>0.001] %>% names()
    ), "white")[1]
    
    
    
    # for perceptual lightness/contrast:
    imhex = imcol %>% hex2RGB() 
    rm(imcol) 
    hexlab = imhex %>% as("LAB")
    lb =  hexlab %>%  coords() %>% .[,1]    # extract L channel
    
    # color range calc (see function for details)
    colorrange = colorfulness(imhex)
    rm(imhex)
    colorrange2 = colorfulness2(hexlab)
    rm(hexlab)
    
    contrastsd =  sd(lb/100) 
    contrastsd=ifelse(is.na(contrastsd), 0, contrastsd)
    contrastrange = table(round(lb)) %>%
      .[.>params$pixeltotal*0.001] %>% # exclude outliers (few super bright or dark pixels)
      names() %>% as.numeric() %>% range() %>% 
      {abs(.[2]-.[1])/100} # standardize to 0-1
    rm(lb)
    gc()

    
    
    
    sizes = vector("list", length(coefs)) #; names(sizes)=as.character(coefs)
    for(i in seq_along(coefs) ){
      if(i==1){      # first element is hard-coded to be baseline size
        x = baseline
        bmpsize  = image_write(x, raw(), "rgb", depth=16) %>% length() 
        baseinisize = bmpsize # will do other compressions in the end which will use this
        tmp = transformer(x, fillcol=fillcol)
      } else { 
        x = image_resize(baseline, paste0(newsize[[i]][1],"x", newsize[[i]][2], "!") )
        tmp = transformer_small(x, fillcol=fillcol) # resized images only transformed with a smaller subset (to save time, and also for half the transformations the size doesn't seem to make any difference)
      }
      
      # smaller = more compressible, simpler;
      gftmp=sizer(x, "gif")
      pngtmp=sizer(x, "png", depth=8)
      sgif = sizer(tmp, "gif") / gftmp
      #spng = sizer(tmp, "png", depth=16) / pngtmp
      
      names(sgif) = paste0(names(sgif), "_gif")
      #names(spng) = paste0(names(spng), "_png")
      
      # compressions without transforms but different compressors
      notr = c(
        compress_gif=gftmp,    # itself, gif; 8bit by default i.e. 256 cols
        compress_png=pngtmp,    
        # bmp = sizer(x, "bmp") ,    # bitmap is smaller than rgb raw
        compress_jpeg0=sizer(x, "jpeg", quality=0),   # jpeg 0 is half the size of 100, lossy
        compress_jpeg100=sizer(x, "jpeg", quality=100) 
        # tiff_lzw=sizer(x, "tiff",depth=8, compression="LZW"), # lzw sometimes crashes for no clear reason (something w image tags), so just leave it out
        #tiff_zip=sizer(x, "tiff",depth=8, compression="ZIP"), # also buggy
        #svg=sizer(x, "svg",depth=16), # same as png, probably just embeds it
      ) / baseinisize   # normalize by bitmap/uncompressed size
      
      sizes[[i]] = c(
        sgif, 
        notr,
        lines_edge5_gray_png = sizer(tmp$lines_edge5_gray, "png", depth=8) / pngtmp,
        blur10_png = sizer(tmp$blur10, "png", depth=8) / pngtmp,
        lines_division_gray_png = sizer(tmp$lines_division_gray, "png", depth=8) / pngtmp
      )
      
      # do fourier transforms - needs to be resized separately; will do only for baseline
      
      if(i==1){
        # only do with baseline size
        forms = c("colors_quantize3", "blur10", "lines_division_gray")
        fftsizes = c(
          fourier_sizer(baseline),
          fourier_sizer(tmp[[forms[1]]]),
          fourier_sizer(tmp[[forms[2]]]),
          fourier_sizer(tmp[[forms[3]]])
        ) / baseinisize
        names(fftsizes) = 
          c("fft1", "fft2",paste0(c("fft1", "fft2"), "_", rep(forms,each=2)))
      }
      
      # notransform compression size instead of bmp, as all imgs have some baseline compressibility to them.
      
      names(sizes[[i]]) = paste0(names(sizes[[i]]), "_x", coefs[i])
      
      # image_animate(tmp %>% image_join(), delay=70)
      # image_animate(tmp[order(sizer(tmp)/inisize, decreasing = T)] %>% image_join(), delay=70)
      # print(paste(Sys.time(),i, "coef done"))
    }
    sizes2 = unlist(sizes, F, use.names = T) # to vector
    # standardized unique names now
    
    
    # do few more extra small conversions with just the compressor, no transforms
    extrasizes=NULL
    if(!is.null(extracoefs)){
      extrasizes = vector("list", length(extracoefs))
      for(i in seq_along(extracoefs)){
        newsize2 =  round(xy * sqrt( (pixeltotal*(extracoefs[i]^2)) / (xy[1] * xy[2])))
        x2 = image_resize(baseline, paste0(newsize2[1],"x", newsize2[2], "!") )
        extrasizes[[i]] = c(
          gif=sizer(x2, "gif"),    # itself, gif; 8bit by default i.e. 256 cols
          png=sizer(x2, "png", depth=8),    # itself
          #tiff_lzw=sizer(x2, "tiff",depth=8, compression="LZW"),
          #tiff_zip=sizer(x2, "tiff",depth=8, compression="ZIP"),
          #svg=sizer(x2, "svg",depth=8),
          jpeg0=sizer(x2, "jpeg", quality=0) # lossy compression too
        ) /  baseinisize
        names(extrasizes[[i]]) = paste0("compress_", names(extrasizes[[i]]), "_x", extracoefs[i])
      }
      extrasizes = unlist(extrasizes, F, use.names = T) # to vector
    }
    
    
    # add summary color distribution stats (might as well)
    colfreq = count_colors(image_quantize(baseline,200, dither=F, colorspace = "lab"), pixelsize, forplot=F) %>% {.[,2]/sum(.[,2])}
    
    # add entropy of estimated composition lines
    angleentropy = image_despeckle(baseline,2)  %>% image_canny() %>% 
      image_hough_txt("40x40+30")  %>% strsplit("#") %>% unlist(F,F) %>% 
      .[-(1:4)] %>% gsub("^ ([0-9]+) .*", "\\1",.) %>% as.numeric() %>% 
      .[!is.na(.)] %>% entropy()  
    # if no lines, returns 0 entropy (noise should yield noisy canny lines > many hough lines > high entropy; but not impossible this fails sometimes, if the canny algo  doesn't pick up any structure at all in the noise)
    
    
    sdcolfreq = sd(colfreq)
    sdcolfreq = ifelse(is.na(sdcolfreq), 0, sdcolfreq)
    compvec = c(
      sizes2        # transforms;  has names
      , extrasizes  # more sizes no compression
      , fftsizes    # fourier on baseline & 3 transforms
      , fdcalc(baseline)  # 3 fractal dim estimates,  have names
      , stats_colfreq_mean = 1-mean(colfreq)
      , stats_colfreq_median = 1-median(colfreq) # mean/med/max: inverse, higher val=less colorful
      , stats_colfreq_sd = 1-sdcolfreq
      , stats_colfreq_max = 1-max(colfreq)
      , stats_colfreq_entropy = entropy::entropy(colfreq)
      , stats_colorfulness_rgb = colorrange
      , stats_colorfulness_lab = colorrange2
      , stats_contrastsd = contrastsd
      , stats_contrastrange = contrastrange
      , stats_angleentropy = angleentropy
    )
    attr(compvec, "fp") = fp  # store filepath for debug, matching, and general sanity check
    
    try({
      # cleanup - otherwise bigger machine will run into mem issues
      image_destroy(x)
      image_destroy(baseline)
      tmp=lapply(tmp, image_destroy)
      image_destroy(x2) # only exists if extracoefficients
      
    })
    rm(list=setdiff(ls(all.names = TRUE), "compvec" ) )
    gc(verbose = F, full=T)
    gc(verbose = F, full=T)
  }, 
  error=function(e){errorprinter(paste(e, fp))},
  warning=function(e){errorprinter(paste(e, fp))} )
  return(compvec)
}


docompressions_bw = function(fp, params){
  # use image_crop or image_border to get center
  compvec=NA # if pipeline fails at any point, will return 1-length NA vector, easy to catch later
  tryCatch({
    # pixeltotal=250000; pixelsize=500
    pixeltotal=params$pixeltotal
    pixelsize =params$pixelsize
    coefs = params$coefs
    extracoefs = params$extracoefs
    #normsize = params$normsize
    #mainalgo = params$mainalgo
    minsize = params$minsize
    # if(is.null(normsize)){normsize="gif"} # if not specified, assume normalizer is baseline size
    # if(is.null(minsize)){minsize=0.5} # if not specified, assume size threshold is half of baseline
    # if(is.null(mainalgo)){mainalgo="gif"}  # by default, gif
    
    im = image_read(fp) %>% image_convert(format="rgb", depth = 8, antialias = F)
    
    im = im %>% image_quantize(2, dither = F, colorspace = "Gray") %>%  image_fx(expression = "round(p)")  # convert to black and white if already wasn't (if has some dark grays in it; would produce unexpected results on full color images); quantize step makes them a bit smoother; resulting colorspace is gray though, not actual bw.
    
    xy = as.numeric(image_info(im)[2:3])
    
    # if input image too small, fail with a message, caught by trycatch; vector will be NA
    if(prod(xy) < (pixeltotal*minsize) ){
      stop(paste0("Image too small ", xy[1],"x", xy[2]) )
    }
    
    
    # calculate best x*y ratio to get equivalent total px for baseline and then downsample
    newsize = list()
    for(i in seq_along(coefs)){
      if(coefs[i]==1){
        newsize[[i]]=round(xy * sqrt(pixeltotal / (xy[1] * xy[2])))
      } else {
        newsize[[i]] =  round(xy * sqrt( (pixeltotal*(coefs[i]^2)) / (xy[1] * xy[2])))
      }
    }
    newsize = newsize %>% lapply(., function(x) if(all(x <= xy)){x}else{xy}   )  # but make sure no upsampling
    
    # print(paste(Sys.time(), "new sizes calculated"))
    
    if(length(newsize)>1){
      # second safeguard: will fail if file too small to be downsampled comparably
      if(all(newsize[[1]] <= newsize[[2]])){
        warning(paste(fp, "too small native resolution"))
        return(NA) # returns single NA
      }
    }
    # should check lengths and set row to be NAs after the parLapply part
    
    if(all(newsize[[1]] == xy ) ){
      # if already exactly correct size (e.g. the norms), then skip the unnecessary resize step
      baseline = im %>%  image_flatten() 
    } else {
      baseline = image_resize(im, paste0(newsize[[1]][1],"x", newsize[[1]][2], "!") ) %>% 
        image_flatten() # flattening by itself to remove transparency
    }
    rm(im)
    
    # get col triplets; will reuse multiple places below
    #imcol = image_raster(baseline)$col
    
    # infer suitable color for flood-filling, so that it doesn't overlap with the image's primary colors; neon green or magenta probably least likely; but include backup white if none fit (might be some very neon/primary color image which would have all these)
    fillcol = rgb(0,1,0)
    
    sizes = vector("list", length(coefs)) #; names(sizes)=as.character(coefs)
    for(i in seq_along(coefs) ){
      if(i==1){      # first element is hard-coded to be baseline size
        x = baseline
        bmpsize  = image_write(x, raw(), "rgb", depth=16) %>% length()
        baseinisize = bmpsize # will do other compressions in the end which will use this
        tmp = transformer_bw(x, fillcol=fillcol)
      } else { 
        stop("bw transformer meant for only one size")
        #x = image_resize(baseline, paste0(newsize[[i]][1],"x", newsize[[i]][2], "!") )
        #tmp = transformer_small(x, fillcol=fillcol) # resized images only transformed with a smaller subset (to save time, and also for half the transformations the size doesn't seem to make any difference)
      }
      
      # smaller = more compressible, simpler;
      gftmp=sizer(x, "gif")
      pngtmp=sizer(x, "png", depth=8)
      sgif = sizer(tmp, "gif") / gftmp
      #spng = sizer(tmp, "png", depth=16) / pngtmp
      
      names(sgif) = paste0(names(sgif), "_gif")
      #names(spng) = paste0(names(spng), "_png")
      
      # compressions without transforms
      notr = c(
        gif=gftmp,    # itself, gif; 8bit by default i.e. 256 cols
        png=pngtmp,    
        bmp = sizer(x, "bmp") ,    # bitmap is smaller than rgb raw
        jpeg0=sizer(x, "jpeg", quality=0),   # jpeg 0 is half the size of 100, lossy
        jpeg100=sizer(x, "jpeg", quality=100) 
        # tiff_lzw=sizer(x, "tiff",depth=8, compression="LZW"), # lzw sometimes crashes for no clear reason (something w image tags), so just leave it out
        #tiff_zip=sizer(x, "tiff",depth=8, compression="ZIP"), # also buggy
        #svg=sizer(x, "svg",depth=16), # same as png, probably just embeds it
      ) / baseinisize   # normalize by bitmap/uncompressed size
      
      
      sizes[[i]] = c(sgif, 
                     notr,
                     lines_edge5_gray_png = sizer(tmp$lines_edge5_gray, "png", depth=8) / pngtmp,
                     blur10_png = sizer(tmp$blur10, "png", depth=8) / pngtmp
                     # lines_division_gray_png = sizer(tmp$lines_division_gray, "png", depth=8) / pngtmp
      )
      
      # notransform compression size instead of bmp, as all imgs have some baseline compressibility to them.
      
      names(sizes[[i]]) = paste0(names(sizes[[i]]), "_x", coefs[i])
      
      # image_animate(tmp %>% image_join(), delay=70)
      # image_animate(tmp[order(sizer(tmp)/inisize, decreasing = T)] %>% image_join(), delay=70)
      # print(paste(Sys.time(),i, "coef done"))
    }
    sizes2 = unlist(sizes, F, use.names = T) # to vector
    # standardized unique names now
    
    
    # do few more extra small conversions with just the compressor
    extrasizes=NULL
    if(!is.null(extracoefs)){
      extrasizes = vector("list", length(extracoefs))
      for(i in seq_along(extracoefs)){
        newsize2 =  round(xy * sqrt( (pixeltotal*(extracoefs[i]^2)) / (xy[1] * xy[2])))
        x2 = image_resize(baseline, paste0(newsize2[1],"x", newsize2[2], "!") )
        extrasizes[[i]] = c(
          gif=sizer(x2, "gif"),    # itself, gif; 8bit by default i.e. 256 cols
          png=sizer(x2, "png", depth=8),    # itself
          #tiff_lzw=sizer(x2, "tiff",depth=8, compression="LZW"),
          #tiff_zip=sizer(x2, "tiff",depth=8, compression="ZIP"),
          #svg=sizer(x2, "svg",depth=8),
          jpeg0=sizer(x2, "jpeg", quality=0) # lossy compression too
        ) /  baseinisize
        names(extrasizes[[i]]) = paste0(names(extrasizes[[i]]), "_x", extracoefs[i])
      }
      extrasizes = unlist(extrasizes, F, use.names = T) # to vector
    }
    
   
    
    # add entropy of estimated composition lines
    angleentropy = image_despeckle(baseline,2)  %>% image_canny() %>% 
      image_hough_txt("40x40+30")  %>% strsplit("#") %>% unlist(F,F) %>% 
      .[-(1:4)] %>% gsub("^ ([0-9]+) .*", "\\1",.) %>% as.numeric() %>% 
      .[!is.na(.)] %>% entropy()  
    # if no lines, returns 0 entropy (noise should yield noisy canny lines > many hough lines > high entropy; but not impossible this fails sometimes, if the canny algo  doesn't pick up any structure at all in the noise)
    
    
    
    compvec = c(
      sizes2        # transforms;  has names
      , extrasizes  # more sizes
      , fdcalc(baseline)  # 3 fractal dim estimates,  have names
      # , colfreq_mean = 1-mean(colfreq)
      # , colfreq_median = 1-median(colfreq) # mean/med/max: inverse, higher val=less colorful
      # , colfreq_sd = 1-sd(colfreq)
      # , colfreq_max = 1-max(colfreq)
      # , colfreq_entropy = entropy::entropy(colfreq)
      # , colorfulness_rgb = colorrange
      # , colorfulness_lab = colorrange2
      # , contrastsd = contrastsd
      # , contrastrange = contrastrange
      , angleentropy = angleentropy
    )
    attr(compvec, "fp") = fp  # store filepath for debug, matching, and general sanity check
    print(paste(Sys.time(), fp))
    
  }, error=function(e){try(print(paste(fp, e)))}, warning=function(e){try(print(paste(fp, e)))} )
  return(compvec)
}


errorprinter=function(e, file=file.path(params$projectfolder, "log2.txt")){
  try({
  print(e)
  cat(paste("\n",Sys.time(), e,"\n"),file=file,append=TRUE)
  })
}



# main function, runs in parallel
docompressions_multicore = function(d, params, pathvar="Path2"){
  bw = params$bw
  if(is.null(bw)){bw=F}
  ntry=0
  # try and try again too in case fails, sometimes does for some reason
  while(!exists("compvecs", envir = environment() ) & ntry<4){
    ntry=ntry+1
    nc = detectCores() - params$nfree
    cl = makeCluster(nc, outfile=file.path(params$projectfolder, "log.txt"))
    print(paste(Sys.time(), nc, "cores, starting"))
    tryCatch({
      clusterExport(cl, "params", envir=environment() )
      clusterExport(cl, c("docompressions", "transformer","transformer_small",
                          "sizer","fourier_sizer", "fdcalc",
                          "count_colors","scramble_colors", "circleFun","colorfulness","colorfulness2",
                          "docompressions_bw", "transformer_bw", "errorprinter"
      ),envir = globalenv() ) 
      
      clusterEvalQ(cl, suppressPackageStartupMessages(c(library("magick"), library("dplyr"), library("entropy"), library("gifski"), library("ggplot2"), library("colorspace"), library("fractaldim"), library("reshape2")) ) )
      clusterEvalQ(cl, magick:::magick_threads(1)) # unexported magick function, if this is not set, then magick will also try to use all cores, which can quickly lead to memory problems on larger machines with 20+ cores.
      
      # do computation:
      compvecs = parLapply(cl, d[[pathvar]],  docompressions, params=params ) 
    }, 
    error=function(e){errorprinter(e)},
    warning=function(e){errorprinter(e)},  
    finally = stopCluster(cl) 
    )
  }
  
  try(stopCluster(cl)) # try just in case
  try(rm(cl))
  try(closeAllConnections()) # in case the stopcluster left something behind or failed
  gc(verbose = F, full=T)
  
  if(!exists("compvecs")){
    errorprinter("---------failure---------")
    return(NULL)
  } else {
    
  }
  tryCatch({
  errorprinter( 
    paste( "done,", sum(sapply(compvecs, length)>1 ), "/", nrow(d), "probably succeeded")
    )
    },error=function(e){errorprinter(e)})
  return(compvecs)
}




