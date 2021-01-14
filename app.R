# Copyright (C) 2021, Phebo Wibbens
  
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(shiny)
library(shinyMatrix)
library(tidyverse)

precision <- 1e-3 
nV <- 2
m <- matrix(0, nrow = nV, ncol = 2)
rownames(m) <- paste("Trial", 1:nV)
colnames(m) <- c("Vaccinated", "Control")

post <- function(alpha, nv, nc) {
  p <- (1-alpha)^nv / (2-alpha)^{nv+nc}
  p / sum(p)
}

ui <- fluidPage(
  titlePanel("Vaccine efficacy"),
  sidebarLayout(
    sidebarPanel(
      
      matrixInput("m", label = "Enter the number of trial participants with a particular outcome (e.g. contracted COVID) in vaccinated and control groups",
                  value = m, class = "numeric", rows = list(editableNames = T, names = T, extend = T),
                  cols = list(names = T), copy = T, paste = T)
    ),
    mainPanel(plotOutput("plot"))
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    df <- as_tibble(input$m) %>% mutate(trial = factor(rownames(input$m), levels=rownames(input$m))) %>%
      select(trial, nv=Vaccinated, nc=Control)
    df2 <- expand_grid(df, alpha = seq(0, 1, precision)) %>% group_by(trial) %>% mutate(
      p = post(alpha, nv, nc),
      cump = cumsum(p),
      pdens = p / precision) %>%
      ungroup()
    df3 <- df2 %>% group_by(trial) %>% summarize(
      median = alpha[which.min(abs(cump-0.5))],
      p.025 = alpha[which.min(abs(cump-0.025))],
      p.975 = alpha[which.min(abs(cump-0.975))])
    ggplot(df3, aes(x = fct_rev(trial), y = median, ymin = p.025, ymax = p.975)) + geom_crossbar() +
      coord_flip() + xlab("") + geom_hline(yintercept = c(0,1), lty=2) +
      scale_y_continuous(labels = scales::percent_format(), name = "Efficacy (median and 95% interval)") +
      theme(text = element_text(size=20))
    })
}

shinyApp(ui = ui, server = server)
