## Color ramp from ggsci but the package wasn't convenient enough to use
## the way I wanted
## A version of simpsonCols with extended color slots for subwatersheds in
## Opunohu and Pao Pao
simpsonCols <- c("Atiha" = "#FED439", 
                 "Haapiti" = "#9b5daa",
                 "Ha'apiti" = "#9b5daa",
                 "Haumi" = "#8A9197", 
                 "Maatea" = "#D2AF81",
                 "Ma'atea" = "#D2AF81",
                 "Maharepa" = "#FD7446",
                 "Opunohu" = "#bdd56c", 
                 "Opunohu 1" = "#D5E4A2", "Opunohu 2" = "#acca47",
                 "Paopao" = "#1E90dd",
                 "Paopao 1" = "#4aa9e7", "Paopao 2" = "#197EC0", "Paopao 3" = "#125c8c", 
                 "Papetoai" = "#F05C3B",
                 "Pihaena" = "#91331F",
                 "Teavaro" = "#46732E", 
                 "Vaianae" = "#71D0F5")

simpsonColsFill = col2rgb(simpsonCols) %>%
  apply(MARGIN = 2,
        FUN = function(x){
          x + 0.5*(255-x)
        }) %>%
  apply(MARGIN = 2,
        FUN = function(x){
          rgb(x[1],x[2],x[3],
              maxColorValue = 255)
        })
