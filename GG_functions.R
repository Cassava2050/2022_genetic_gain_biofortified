library(lmerTest)

# -------------------------------------------------------------------------
# utilities ---------------------------------------------------------------
# -------------------------------------------------------------------------

library(lucid)

Heri.asreml <- function(Model,genotype){
  
  # Genetic variance component
  vc.g <- vc(Model)[grepl(genotype,vc(Model)$effect),2]
  vc.g 
  
  # Mean variance of a difference of two genotypic BLUPs
  vdBLUP.mat <- predict(Model, classify=genotype, only=genotype, sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
  vdBLUP.avg 
  
  #############
  # H2 Cullis #
  #############
  H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
  
  return(H2Cullis)
  
}

Heri.AR.AR <- function(Model,genotype){
  sp4.pv.h2 <- predict(Model, classify = genotype, only = genotype, sed = TRUE)
  avepev <- mean(sp4.pv.h2$pvals$std.error^2)
  vcV <- vc(Model)[grepl(genotype,vc(Model)$effect),2]
  h2.wn4 <- 1 - avepev/(vcV); h2.wn4[1]
  h2.wn4
}


# -------------------------------------------------------------------------



ran <- function(var){
  effect <- paste0("(",1,"|",var,")")
  return(effect)
}

VarG <- function(model, comp){
  v <- as.data.frame(VarCorr(model))
  v <- v[v$grp==comp,"vcov"]
  return(v)
}

VarE <- function(model){
  v <- as.data.frame(VarCorr(model))
  v <- v[v$grp=="Residual","vcov"]
  return(v)
}

h.cullis <- function(model, gen){
  aveped <- mean(attr(ranef(model,drop=T)[[gen]],"postVar"))
  vc.g <- as.data.frame(VarCorr(model))
  vc.g <- vc.g[vc.g$grp==gen,"vcov"]
  ifelse(vc.g==0, 0 , round(1-aveped/vc.g,3) )
}

varG.pvalue <- function(model, gen){
  table <- try(suppressWarnings(broom.mixed::tidy(lmerTest::ranova(model))), silent = T)
  if(length(class(table))==1){
    return(NA)
  } else{
    term <- grepl(gen, x = table$term)
    as.numeric(table[term,"p.value"])
  }
}

lme4_res <- function(model, return=F){
  res <-  residuals(model, scaled=TRUE)
  data <- model@frame
  data$res <- res
  data$Classify <- NA
  data$Classify[which(abs(data$res)>=3)] <- "Outlier" 
  data$Classify[which(abs(data$res)<3)  ] <- "Normal"
  ix = ifelse(length(which( abs(res)>3 ))>=1, length(which( abs(res)>3 )) , 0  )
  if (return) {
    return(data)
  } else{
    return(ix)
  }
}

mult_summary <- function(models, gen = "Name", y = "response"){
  exp <- names(models)
  gv <- unlist(lapply(models, VarG, gen))
  ev <- unlist(lapply(models, VarE))
  he <- unlist(lapply(models, h.cullis, gen ))
  out <- unlist(lapply(models, lme4_res ))
  summ <- data.frame(Experiment=exp, y = y ,varG = gv, varE = ev, h2 = he, outliers=out , row.names = NULL)
  return(summ)
}

mult_lme4 <- function(data, equation, var_sub ){
  models <- list()
  data[,var_sub] <- as.factor(data[,var_sub])
  for (exp in levels(data[,var_sub])) {
    tmp_dt <- dplyr::filter(data,.data[[var_sub]]%in%exp)
    model <-  try(lmerTest::lmer(equation,data=tmp_dt, na.action = na.omit), silent = T)
    if (class(model)=="try-error") {
      models[[exp]] <- NULL
    } else{
      models[[exp]] <- model
    }
  }
  return(models)
}


lme4_BLUPs <- function(model, genotype){
  BLUPS <- ranef(model)[[genotype]]
  BLUPS <- data.frame(as.factor(row.names(BLUPS)),BLUPS[,1])
  colnames(BLUPS) <- c("Genotype","Effect")
  BLUPS <- dplyr::arrange(BLUPS,desc(Effect))
  BLUPS <- data.frame(BLUPS[,1],round(BLUPS[,2],2))
  names(BLUPS) <- c("Line","BLUPs")
  d <- broom.mixed::augment(ranef(model))
  d <- d[d$grp==genotype,c("level","std.error")]
  d <- data.frame(level=d[,1],std.error=round(d[,2],2))
  BLUPS <- merge(BLUPS,d,by.x="Line",by.y="level")
  BLUPS
}


lme4_plotly <- function(blups){
  BLUPS <- blups
  BLUPS$Lu <- BLUPS[,2]-1.645*BLUPS[,3]
  BLUPS$Ls <- BLUPS[,2]+1.645*BLUPS[,3]
  v <- as.character(BLUPS[order(BLUPS[,2],decreasing = TRUE),1])
  names(BLUPS)[2] <- "predicted.value"
  p <- ggplot(BLUPS,aes(x=Line , y=predicted.value))+
    geom_point(size = 1) +
    geom_errorbar(aes(ymax = Ls, ymin = Lu))+
    theme_bw() +
    geom_hline(yintercept = mean(BLUPS[,2]), linetype=2 ,color="red")+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    ylab(names(BLUPS)[2])+scale_x_discrete(limits=v)
  plotly::ggplotly(p)
}


lme4_ggplot <- function(blups, title = NULL , subtitle = NULL){
  BLUPS <- blups
  BLUPS$Lu <- BLUPS[,2]-1.645*BLUPS[,3]
  BLUPS$Ls <- BLUPS[,2]+1.645*BLUPS[,3]
  v <- as.character(BLUPS[order(BLUPS[,2],decreasing = TRUE),1])
  names(BLUPS)[2] <- "predicted.value"
  p <- ggplot(BLUPS,aes(x=Line , y=predicted.value))+
    geom_point(size = 1) +
    geom_errorbar(aes(ymax = Ls, ymin = Lu))+
    theme_bw() +
    geom_hline(yintercept = mean(BLUPS[,2]), linetype=2 ,color="red")+
    theme_ipsum(base_size = 10) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1),
          axis.ticks.x=element_blank())+
    labs(x="", y=names(BLUPS)[2],
         title=title,
         subtitle=subtitle) +
    scale_x_discrete(limits=v)
  return(p)
}


mult_asreml <- function(data, genotype, trait ,fixed_comp,  random_comp, comp_equation ,var_sub, clean_out, iterations = 1 ){
  models <- list()
  pvals <- list()
  data <- as.data.frame(data)
  data[,var_sub] <- as.factor(data[,var_sub])
  for (exp in levels(data[,var_sub])) {
    tmp_dt <- dplyr::filter(data,.data[[var_sub]]%in%exp)
    
    equation <- reformulate(comp_equation, response = trait)
    mod0 <-  try(lmerTest::lmer(equation, data = tmp_dt, na.action = na.omit), silent = T)
    
    fixed <- reformulate(fixed_comp, response = trait)
    random <- reformulate(random_comp, response = NULL)
    model <-  try(asreml::asreml(fixed = fixed, random = random ,data=tmp_dt, na.action=list(x="include",y="include"), trace =0), silent = T)
    if (class(model)=="try-error") {
      models[[exp]] <- NULL
    } else{
      models[[exp]] <- model
      
      if(clean_out){
        resum_out <- lme4_res(mod0, return = F)
        if(resum_out>0){
          tmp_out <- 1
          counter <- 1
          while (tmp_out>0 & counter<=iterations) {
            c_datos <- lme4_res(mod0, return = T)
            c_datos[c_datos$Classify=="Outlier", trait] <- NA
            mod0 <-  try(lmerTest::lmer(equation, data = c_datos, na.action = na.omit), silent = T)
            model <-  try(asreml::asreml(fixed, random = random, data = c_datos, na.action=list(x="include",y="include"), trace =0  ), silent = T)
            tmp_out <- lme4_res(mod0, return = F)
            if(iterations>1) resum_out <-  resum_out + tmp_out
            counter <- counter + 1
          }
        }
        models[[exp]] <- model
      }
      
      pvals[[exp]] <- predict(model, classify = genotype)$pvals
    }
  }
  return(list(models=models, pvals = pvals))
}


res_data_lme4 <- function(Model){
  if(class(Model)=="lm"){
    Data <- Model$model
    VarE <- sigma(Model)^2
  } else {
    Data <- Model@frame
    VarE <- VarE(Model) 
  }
  Data$Index <- 1:length(residuals(Model))
  Data$Residuals <- residuals(Model)
  u <- +3*sqrt(VarE)
  l <- -3*sqrt(VarE)
  Data$Classify <- NA
  Data$Classify[which(abs(Data$Residuals)>=u)] <- "Outlier" 
  Data$Classify[which(abs(Data$Residuals)<u)  ] <- "Normal"
  Data$l <- l
  Data$u <- u
  Data$fit <-  fitted.values(Model)
  return(Data)
}

res_qqplot <- function(data_out, title = NULL){
  q <- dplyr::filter(data_out,!is.na(Classify)) %>% 
    ggpubr::ggqqplot(x="Residuals",
                     fill="Classify",
                     ggtheme=theme_ipsum(),
                     ylab = "Sample Quantile", 
                     xlab = "Theoretical Quantile", title =title )
  return(q)
}


# MET ---------------------------------------------------------------------

var_fa <- function(model){
  ASM <- MrBean:::fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  psi = ASM$gammas[[1]]$`specific var`
  VarTot = sum(diag(L.star %*% t(L.star))) / sum(diag(L.star %*% t(L.star) + diag(psi) ))
  
  paf.site <- ASM$gammas[[1]]$`site %vaf`
  
  VarGenEnv <- diag(L.star %*% t(L.star) + diag(psi) )
  TotGenVar <- sum(VarGenEnv)
  
  VarFA1 <- sum(VarGenEnv*paf.site[,1])/100
  VarFA2 <- sum(VarGenEnv*paf.site[,2])/100
  
  PerVarFA1 <- round(VarFA1/TotGenVar*100,1)
  PerVarFA2 <- round(VarFA2/TotGenVar*100,1)
  
  return(perc = c(PerVarFA1,PerVarFA2))
}


"circleFun" <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

require(ggrepel)

biplot_fa2 <- function(model,
                       predictions,
                       one = -1, 
                       second = -1,
                       fscore = 2, 
                       fmult = 10,
                       alpha = 0.3,
                       alpha_site = 0.5, 
                       alpha_ind = 0.5,
                       subtitle=NULL, 
                       gen = NULL,
                       size_ind_biplot = 3){
  
  ASM <- MrBean:::fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  L.star[,1] <- L.star[,1]*one
  L.star[,2] <- L.star[,2]*second
  psi = ASM$gammas[[1]]$`specific var`
  Gvar <- ASM$gammas[[1]]$Gmat
  Cmat <- ASM$gammas[[1]]$Cmat
  Env.means <- predictions
  names(Env.means)[1:2] <- c("site", "BLUE") 
  faComp <- data.frame(site = rownames(L.star), fa1 = L.star[,1], fa2 = L.star[,2], psi = psi, Vg = diag(Gvar), BLUE = Env.means$BLUE)
  
  percentg <- var_fa(model)
  
  # Without Standardize Loadings
  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  loadings = ggplot(faComp, aes(x = fa1, y = fa2)) + 
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") + 
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) + 
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha) 
  
  # Standardize Loadings
  faCompR <- faComp
  faCompR[,2:3] <- diag(1/sqrt(diag(Gvar))) %*% L.star 
  d <- data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=faCompR[,2], vy=faCompR[,3]) 
  circle <- circleFun(c(0,0),2,npoints = 100)
  loading_C <- ggplot(faCompR, aes(x = fa1, y = fa2)) + 
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") + 
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) + 
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) + 
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) + 
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha=alpha) +
    geom_path(data=circle, aes(x,y))
  
  
  # Biplot
  fa12_scores = ASM$blups[[1]]$scores
  names(fa12_scores)[2:3] <- c("comp", "Genotype")
  fa12_scores %<>% select(-blup) %>%  spread(. ,"comp", value = "blupr")
  names(fa12_scores) = c("Genotype", "fa1", "fa2")
  fa12_scores$fa1 <- fa12_scores$fa1*one
  fa12_scores$fa2 <- fa12_scores$fa2*second
  fa12_scores$Score <-  ifelse(sqrt(fa12_scores$fa1^2+fa12_scores$fa2^2)>fscore,1,0)
  
  if(!is.null(gen)){
    
    fa12_scores[fa12_scores$Genotype %in% gen, "Score" ] <- 1
    
  }
  
  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  biplot = ggplot(faComp, aes(x = fa1, y = fa2)) + 
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") + 
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1,  alpha = alpha_site) + 
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) + 
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha) +
    geom_label_repel(data = subset(fa12_scores, Score==1), 
                     aes(label = Genotype, x = fmult*fa1 , y= fmult*fa2),
                     colour = "red",segment.colour = "red" , size=size_ind_biplot, alpha = alpha_ind) 
  
  
  return(list(g1=loadings, g2= loading_C, g3 = biplot))
}


# ggcor


ggCor <- function(myData, colours = c('#db4437','white','#FF9D00'),
                  blackLabs = c(-0.7, 0.7), showSignif = TRUE,
                  pBreaks = c(0, .0001, .001, .01, Inf), pLabels = c('***','**','*', 'ns'),
                  showDiagonal = FALSE, Diag = NULL, returnTable = FALSE, returnN = FALSE, text_diag = 3){
  
  #   Goal      : Return a ggplot object to plot a triangular correlation figure between 2 or more variables.
  #               Depends on the packages 'ggplot2' 'psych' and 'reshape'
  #
  #   Input     : myData       = A data.frame with numerical columns for each variable to be compared.
  #   Input     : colours      = A vector of size three with the colors to be used for values -1, 0 and 1.
  #   Input     : blackLabs    = A numeric vector of size two, with min and max correlation coefficient 
  #                              limits to display with black tags. Any value outside this range will be 
  #                              displayed with white tags.
  #   Input     : showSignif   = Logical scalar. Display significance values ?
  #   Input     : pBreaks      = Passed to function 'cut'. Either a numeric vector of two or more unique 
  #                              cut points or a single number (greater than or equal to 2) giving the
  #                              number of intervals into which x is to be cut.
  #   Input     : pLabels      = Passed to function 'cut'. labels for the levels of the resulting category.
  #                              By default, labels are constructed using "(a,b]" interval notation. 
  #                              If pLabels = FALSE, simple integer codes are returned instead of a factor.
  #   Input     : showDiagonal = Logical scalar. Display main diagonal values ?
  #   Input     : Diag         = A named vector of labels to display in the main diagonal. The names are 
  #                              used to place each value in the corresponding coordinates of the diagonal.
  #                              Hence, these names must be the same as the colnames of myData
  #   Input     : returnTable  = Return the table to display instead of a ggplot object
  #
  #   Output    : A ggplot object containing a triangular correlation figure with all numeric variables 
  #               in myData. If returnTable is TRUE, the table used to produce the figure is returned instead.
  #   Authors   : darizasu
  #    Last update: Oct 14, 2020
  
  require(ggplot2)
  
  # Drop non numeric columns in the dataset
  if (sum( !sapply(myData, is.numeric) )){
    
    message('Dropping non-numeric columns in the dataset:\n',
            paste(names( which(!sapply(myData, is.numeric)) ),
                  collapse = '\t'))
    
    myData = myData[,sapply(myData, is.numeric)]
  }
  
  # Calculate corr-coeffs and p values
  cors = psych::corr.test(myData, use = 'pairwise.complete.obs')
  
  # Use the adjusted p values for multiple testing instead of raw coeffs
  cors$p = t(cors$p)
  
  # Keep only the matrices with correlation coefficients, p values and N shared samples
  cors = cors[c(1,2,4)]
  
  # Make sure you have a full matrix of N shared samples
  if (is.vector(cors$n))
    cors$n <- matrix(cors$n, ncol(cors$p), nrow(cors$p),
                     dimnames = dimnames(cors$p))
  
  # For each matrix, do ...
  cors = lapply(cors, function(x){
    
    # Keep the upper triangle of the matrix
    x[upper.tri(x)] = NA
    
    # Transpose the matrix to plot the lower triangle
    x  = as.data.frame(t(x))
    
    # Reshape the matrix to tidy format
    x[,'col'] = colnames(x)
    x  = reshape::melt(x, id='col')
    colnames(x) = c('col','row','value')
    
    # Round coefficients
    x$name = round(x$value,2)
    
    # Sort the x axis according to myData column order
    x$col = factor(x$col, levels = colnames(myData))
    
    # Reverse the y axis for a triangle plot from top-left to bottom-right
    x$row = factor(x$row, levels = rev(colnames(myData)))
    
    # Remove NAs
    x = na.omit(x)
    
  })
  
  # Combine both dataframes with p values and corr coefficients
  cors = merge(x = merge(x = cors$r, y = cors$p, by = c('col','row')),
               y = cors$n, by = c('col','row'))
  
  # Keep x, y, p val and corr-coefficients columns
  cors = cors[,c(1,2,4,5,7)]
  
  if (returnN){
    
    if (returnTable) return(cors)
    
    cors$cols = scale(cors$value, center = T, scale = T)
    cors$cols = ifelse(abs(cors$cols) < 2, 'black', 'white')
    
    p = ggplot(data = cors, aes(x = col, y = row, fill = value)) +
      geom_tile(color = 'gray') + labs(x = NULL, y = NULL) + 
      theme_minimal(base_size = 13) +
      geom_text(aes(x = col, y = row, label = value), color = cors$cols, size = 3) +
      scale_fill_gradient(low = colours[2], high = colours[3],
                          limits = c(0, max(cors$value))) + 
      theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = 'none',
            panel.grid.minor.x = element_blank(), panel.grid.major = element_blank())
    
    return(p)
    
  }
  
  if (showSignif){
    
    # Create a categorical variable for p values as defined by pBreaks
    cors$signi = cut(x = cors$value.y,  right = F,
                     breaks = pBreaks, labels = pLabels)
    
    # Join corr-coeff and p-value to display it as a label for each tile
    cors$label = paste(cors$name.x, cors$sign, sep='\n')
    
  } else {
    
    # The label for each tile is the corr-coeff only
    cors$label = cors$name.x
  }
  
  # If there are user-specified values to display in the diagonal
  if (! is.null(Diag)){
    
    # Check the names in Diag are the same than colnames of myData
    if ( sum(! names(Diag) %in% colnames(myData)) ){
      warning("These elements in 'Diag' do not correspond to column names in 'myData':\n",
              paste(names(Diag)[!names(Diag) %in% colnames(myData)],
                    collapse = '\t'))
    }
    
    # The tiles of the diagonal are gray
    cors[cors$col == cors$row, 'name.x'] = NA
    
    # Get the name of x and y levels
    d = as.character(cors[cors$col == cors$row, 'row'])
    
    # Modify the elements of the diagonal and make sure they are displayed
    cors[cors$col == cors$row, 'label'] = Diag[d]
    showDiagonal = TRUE
  }
  
  # Remove the elements of the main diagonal if you don't want to display
  if (!showDiagonal)  cors = cors[cors$col != cors$row,]
  
  # Show darker tiles with white labels for clarity
  cors$txtCol = ifelse(cors$name.x > blackLabs[1] & 
                         cors$name.x < blackLabs[2], 'black', 'white')
  
  # Do not show tile labels for empty tiles.
  # Make tile labels of the diagonal white
  cors$txtCol[is.na(cors$txtCol)] = 'white'
  
  if (returnTable) return(cors)
  
  
  p = ggplot(data = cors, aes(x = col, y = row, fill = name.x)) + 
    geom_tile(color = 'gray') + labs(x = NULL, y = NULL) + theme_minimal(base_size = 13) +
    geom_text(aes(x = col, y = row, label = label), color = cors$txtCol, size = text_diag) +
    scale_fill_gradient2(low = colours[1], mid = colours[2], high = colours[3]) + 
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = 'none',
          panel.grid.minor.x = element_blank(), panel.grid.major = element_blank())
  
  return(p)
}


# model <- m_vc_rr2
# gen = "genotype"
# env = "trial"
# vc.model = "diag"

extractG <- function(model, gen = "genotype" , env = "trial", vc.model = "corv"){
  
  sites <- data.frame(model$mf)[, env]
  s <- nlevels(sites)
  
  vc <- summary(model)$varcomp
  VCOV <- matrix(0, ncol = s, nrow = s)
  CORR <- matrix(0, ncol = s, nrow = s)
  diag(CORR) <- rep(1, s)
  
  gxe <- paste(env, gen, sep = ":")
  
  if (vc.model == "diag") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    diag(VCOV) <- vc[, 1]
  }
  if (vc.model == "corv") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- rep(vc[2, 1], s)
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "corh") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- vc[2:(s + 1), 1]
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "corgv") {
    vc.corr <- vc[grep(".cor", rownames(vc)), ]
    vc.var <- vc[-grep(".cor", rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        if (i != j) {
          CORR[i, j] <- vc.corr[k, 1]
          CORR[j, i] <- vc.corr[k, 1]
          k <- k + 1
        }
      }
    }
    D <- rep(vc.var[1, 1], s)
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "fa1") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    R <- vc.var[, 1]
    L <- vc.fa1[, 1]
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "fa2") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L <- cbind(L1, L2)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "fa3") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    vc.fa3 <- vc[grep("!fa3", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L3 <- vc.fa3[, 1]
    L <- cbind(L1, L2, L3)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "fa4") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    vc.fa3 <- vc[grep("!fa3", rownames(vc)), ]
    vc.fa4 <- vc[grep("!fa4", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L3 <- vc.fa3[, 1]
    L4 <- vc.fa4[, 1]
    L <- cbind(L1, L2, L3, L4)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "corgh") {
    vc.corr <- vc[grep(".cor", rownames(vc)), ]
    vc.var <- vc[-grep(".cor", rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        if (i != j) {
          CORR[i, j] <- vc.corr[k, 1]
          CORR[j, i] <- vc.corr[k, 1]
          k <- k + 1
        }
      }
    }
    D <- vc.var[1:s, 1]
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "us") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        VCOV[i,j] <- vc[k,1]
        k <- k+1
      }
    }
    VCOV[upper.tri(VCOV)] = t(VCOV)[upper.tri(VCOV)]
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "rr2") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L <- cbind(L1, L2)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  colnames(VCOV) <- levels(sites)
  colnames(CORR) <- levels(sites)
  rownames(VCOV) <- levels(sites)
  rownames(CORR) <- levels(sites)
  
  return(list(VCOV = VCOV , CORR = CORR, vc.model = vc.model))
  
}
