# Weight distribution evolution per census year
p<-ggplot(SHEEP,aes(exp(size)))+geom_histogram(aes(y=..density..),fill="blue",alpha=0.2)+geom_density() + facet_grid(sheep.yr~.)
ggsave("TimeSeriesDensityHist.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 8,height = 5)
# Growth:
p<-ggplot(SHEEP,aes(exp(prev.size),exp(size)))+stat_density2d_filled(contour_var = "ndensity") + geom_abline(slope=1) + xlim(c(10,30)) + ylim(c(10,30)) +
  labs(x="Previous Weight [kg]",y="Weight [kg]",fill="Normalised Density",title="Growth Kernel") +theme(plot.title = element_text(hjust = 0.5))
ggsave("GrowthDensity2D.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 4)
# Size of Offspring
p<-ggplot(SHEEP,aes(exp(size),exp(rec1.wt)))+stat_density2d_filled(contour_var = "ndensity") + geom_abline(slope=1) + xlim(c(10,30)) + ylim(c(5,20)) +
  labs(x="Parent Weight [kg]",y="Offspring Weight [kg]",fill="Normalised Density",title="Growth Kernel") +theme(plot.title = element_text(hjust = 0.5))
ggsave("OffspringSizeDensity2D.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 4)
# Number of Offspring
SHEEP$num.born<-SHEEP$off.born ; SHEEP$num.born[SHEEP$num.born<1] <- NA ; SHEEP$num.born%<>%factor()
p<-ggplot(filter(SHEEP,!is.na(num.born)),aes(num.born))+geom_bar(fill="darkorchid3",na.rm = T) + labs(x="Number of Offspring",y="Count")
ggsave("NumOffspring.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 5,height = 4)
# Survival Probability
p<-ggplot(SHEEP,aes(factor(survived),exp(prev.size)))+geom_violin(aes(fill=factor(survived))) + 
  labs(x="Survived",y="Size at Previous Census [kg]")+ scale_fill_discrete(name = "Survived")
ggsave("SurvivalViolin.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 4)
# Reproduction Probability
p<-ggplot(SHEEP,aes(factor(reproduced),exp(size)))+geom_violin(aes(fill=factor(reproduced))) + 
  labs(x="Reproduced",y="Size at Current Census [kg]") + scale_fill_discrete(name = "Reproduced")
ggsave("ReproducedViolin.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 4)

#%%% Lexis diagrams of the population %%%#
devtools::install_github("ottlngr/LexisPlotR")
library(LexisPlotR)

plot_discrete_cbar = function(
  breaks, # Vector of breaks. If +-Inf are used, triangles will be added to the sides of the color bar
  palette = "Greys", # RColorBrewer palette to use
  colors = RColorBrewer::brewer.pal(length(breaks) - 1, palette), # Alternatively, manually set colors
  direction = 1, # Flip colors? Can be 1 or -1
  spacing = "natural", # Spacing between labels. Can be "natural" or "constant"
  border_color = NA, # NA = no border color
  legend_title = NULL,
  legend_direction = "horizontal", # Can be "horizontal" or "vertical"
  font_size = 5,
  expand_size = 1, # Controls spacing around legend plot
  spacing_scaling = 1, # Multiplicative factor for label and legend title spacing
  width = 0.2, # Thickness of color bar
  triangle_size = 0.000001 # Relative width of +-Inf triangles
) {
  require(ggplot2)
  if (!(spacing %in% c("natural", "constant"))) stop("spacing must be either 'natural' or 'constant'")
  if (!(direction %in% c(1, -1))) stop("direction must be either 1 or -1")
  if (!(legend_direction %in% c("horizontal", "vertical"))) stop("legend_direction must be either 'horizontal' or 'vertical'")
  breaks = as.numeric(breaks)
  new_breaks = sort(unique(breaks))
  if (any(new_breaks != breaks)) warning("Wrong order or duplicated breaks")
  breaks = new_breaks
  if (class(colors) == "function") colors = colors(length(breaks) - 1)
  if (length(colors) != length(breaks) - 1) stop("Number of colors (", length(colors), ") must be equal to number of breaks (", length(breaks), ") minus 1")
  if (!missing(colors)) warning("Ignoring RColorBrewer palette '", palette, "', since colors were passed manually")
  
  if (direction == -1) colors = rev(colors)
  
  inf_breaks = which(is.infinite(breaks))
  if (length(inf_breaks) != 0) breaks = breaks[-inf_breaks]
  plotcolors = colors
  
  n_breaks = length(breaks)
  
  labels = breaks
  
  if (spacing == "constant") {
    breaks = 1:n_breaks
  }
  
  r_breaks = range(breaks)
  
  cbar_df = data.frame(stringsAsFactors = FALSE,
                       y = breaks,
                       yend = c(breaks[-1], NA),
                       color = as.character(1:n_breaks)
  )[-n_breaks,]
  
  xmin = 1 - width/2
  xmax = 1 + width/2
  
  cbar_plot = ggplot(cbar_df, aes(xmin=xmin, xmax = xmax, ymin = y, ymax = yend, fill = factor(color, levels = 1:length(colors)))) +
    geom_rect(show.legend = FALSE,
              color=border_color)
  
  if (any(inf_breaks == 1)) { # Add < arrow for -Inf
    firstv = breaks[1]
    polystart = data.frame(
      x = c(xmin, xmax, 1),
      y = c(rep(firstv, 2), firstv - diff(r_breaks) * triangle_size)
    )
    plotcolors = plotcolors[-1]
    cbar_plot = cbar_plot +
      geom_polygon(data=polystart, aes(x=x, y=y),
                   show.legend = FALSE,
                   inherit.aes = FALSE,
                   fill = colors[1],
                   color=border_color)
  }
  if (any(inf_breaks > 1)) { # Add > arrow for +Inf
    lastv = breaks[n_breaks]
    polyend = data.frame(
      x = c(xmin, xmax, 1),
      y = c(rep(lastv, 2), lastv + diff(r_breaks) * triangle_size)
    )
    plotcolors = plotcolors[-length(plotcolors)]
    cbar_plot = cbar_plot +
      geom_polygon(data=polyend, aes(x=x, y=y),
                   show.legend = FALSE,
                   inherit.aes = FALSE,
                   fill = colors[length(colors)],
                   color=border_color)
  }
  
  if (legend_direction == "horizontal") { #horizontal legend
    mul = 1
    x = xmin
    xend = xmax
    cbar_plot = cbar_plot + coord_flip()
    angle = 0
    legend_position = xmax + 0.1 * spacing_scaling
  } else { # vertical legend
    mul = -1
    x = xmax
    xend = xmin
    angle = -90
    legend_position = xmax + 0.2 * spacing_scaling
  }
  
  cbar_plot = cbar_plot +
    geom_segment(data=data.frame(y = breaks, yend = breaks),
                 aes(y=y, yend=yend),
                 x = x - 0.05 * mul * spacing_scaling, xend = xend,
                 inherit.aes = FALSE) +
    annotate(geom = 'text', x = x - 0.1 * mul * spacing_scaling, y = breaks,
             label = labels,
             size = font_size) +
    scale_x_continuous(expand = c(expand_size,expand_size)) +
    scale_fill_manual(values=plotcolors) +
    theme_void()
  
  if (!is.null(legend_title)) { # Add legend title
    cbar_plot = cbar_plot +
      annotate(geom = 'text', x = legend_position, y = mean(r_breaks),
               label = legend_title,
               angle = angle,
               size = font_size)
  }
  
  cbar_plot
}
# Grid
mylexis <- lexis_grid(year_start = min(SHEEP$sheep.yr,na.rm = T), year_end = max(SHEEP$sheep.yr,na.rm = T)+2, age_start = 1, age_end = 15)
# Function to prepare the data to be in the correct format to plot later
PrepLexisData<-function(input){
  output<-data.frame()
  for(i in 1:nrow(input)){
    
    tmp<-input[i,]
    # Bottom-left
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    # Bottom-right
    tmp$sheep.yr<-tmp$sheep.yr+1
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    # Top
    tmp$age<-tmp$age+1
    output%<>%rbind(cbind(tmp,data.frame(id=i)))
    
    # Then the second triangle component
    # Top-left
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    # Top-right
    tmp$sheep.yr<-tmp$sheep.yr+1
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    # Bottom-left
    tmp$age<-tmp$age-1
    tmp$sheep.yr<-tmp$sheep.yr-1
    output%<>%rbind(cbind(tmp,data.frame(id=-i)))
    
  }
  # Convert from year to an actual date (it doesn't matter if it is 1st Jan or any other day)
  output$date<-paste0(as.character(output$sheep.yr),"-01-01")
  
  # Create a colour scheme for the number of births and deaths
  n<-20
  tmp <- hist(log(output$total+1),breaks = n,plot = F) ; tmp<-findInterval(log(output$total+1), tmp$breaks)
  cols<-RColorBrewer::brewer.pal(n = 11, name = "Spectral") ; cols<-colorRampPalette(cols)(n)
  output$cols<-cols[tmp]
  
  return(output)
}
# Calculate the number of sheep per age group, per year
alivedata<-SHEEP%>%group_by(age,sheep.yr)%>%summarise(total=as.integer(sum(survived,na.rm = T)))%>%drop_na()%>%as.data.frame()%>%PrepLexisData
# Calculate number of deaths per age group, per year
survdata<-SHEEP%>%group_by(age,sheep.yr)%>%summarise(total=as.integer(sum(abs(survived-1),na.rm = T)))%>%drop_na()%>%as.data.frame()%>%PrepLexisData
# Calculate number of births occurred, per age group, per year
reprdata<-SHEEP%>%group_by(age,sheep.yr)%>%summarise(total=as.integer(sum(abs(off.born-1),na.rm = T)))%>%drop_na()%>%as.data.frame()%>%PrepLexisData
# Fill up the lexis grid!
# Alive
p<-lexis_polygon(lg = mylexis, x = alivedata$date, y = alivedata$age, group = alivedata$id,fill = alivedata$cols) +
  labs(x="Year",y="Age",title = "Number of Living Sheep (log)") +theme(plot.title = element_text(hjust = 0.5))
ggsave("AliveLexis.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6) 
p<-plot_discrete_cbar((unique(as.integer(exp((hist(log(alivedata$total+1),breaks = 10,plot = F))$breaks)))-1), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("AliveLexisColorbar.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6) 
# Survival
p<-lexis_polygon(lg = mylexis, x = survdata$date, y = survdata$age, group = survdata$id,fill = survdata$cols) +
  labs(x="Year",y="Age",title = "Number of Deaths (log)") +theme(plot.title = element_text(hjust = 0.5))
ggsave("SurvLexis.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6) 
p<-plot_discrete_cbar((unique(as.integer(exp((hist(log(survdata$total+1),breaks = 10,plot = F))$breaks)))-1), 
                   spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("SurvLexisColorbar.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6) 
# Reproduction
p<-lexis_polygon(lg = mylexis, x = reprdata$date, y = reprdata$age, group = reprdata$id,fill = reprdata$cols) + labs(x="Year",y="Age",title = "Number of Births (log)") +theme(plot.title = element_text(hjust = 0.5))
ggsave("ReprodLexis.png", plot=p,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6)
q<-plot_discrete_cbar((unique(as.integer(exp((hist(log(reprdata$total+1),breaks = 10,plot = F))$breaks)))-1), 
                      spacing = "constant", palette="Spectral",legend_direction = "vertical",direction = -1)
ggsave("ReprLexisColorbar.png", plot=q,path = paste0(directory,'Plots/Hamish'),width = 6,height = 6) 
