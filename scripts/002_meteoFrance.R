library(tidyverse)
library(rcartocolor)

#### Import and wrangle precip data ----
precip = read_csv(
  "data/MooreaPrecipitation.csv") %>%
  filter(Date > as.Date("2017-12-31"),
         Date < as.Date("2020-01-01")) %>%
  mutate(Shore = case_when(StationName == "AFAREAITU 2" ~ "East",
                           StationName == "HAAPITI 3" ~ "West",
                           StationName == "HAAPITI 5" ~ "West",
                           StationName == "PAOPAO 1" ~ "North",
                           StationName == "PAPETOAI 3" ~ "North"))


precipMonthly = precip  %>%
  group_by(Shore, Date) %>%
  summarize(Precip_mm = mean(Precip_mm)) %>%
  ungroup() %>%
  mutate(year = lubridate::year(Date),
         month = lubridate::month(Date)) %>%
  group_by(Shore, year, month) %>%
  summarize(Precip_mm = sum(Precip_mm)) %>%
  ungroup() %>%
  mutate(Date = as.Date(paste0(year,"-",month,"-01"))) %>%
  select(-year,-month) %>%
  pivot_wider(names_from = Date, values_from = Precip_mm) %>%
  pivot_longer(-Shore, names_to = "Date", values_to = "Precip_mm") %>%
  mutate(Date = as.Date(Date))
write_csv(precipMonthly,"data/monthly_precip.csv")

figure2 = precipMonthly %>%
  ggplot() +
  geom_rect(data = data.frame(xmin = as.Date(c("2017-12-17",
                                               "2019-01-22")),
                              xmax = as.Date(c("2018-03-16",
                                               "2019-03-16")),
                              ymin = c(0,0), 
                              ymax = c(900,900),
                              group = c(1,2)),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = group),
            fill = "#2297e6",
            alpha = 0.5) +
  geom_rect(data = data.frame(xmin = as.Date(c("2018-07-17",
                                               "2019-07-17")),
                              xmax = as.Date(c("2018-09-17",
                                               "2019-09-17")),
                              ymin = c(0,0), 
                              ymax = c(900,900),
                              group = c(1,2)),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = group),
            fill = "#d4b895",
            alpha = 0.5) +
  geom_bar(aes(x = Date,
               y = Precip_mm,
               group = Shore,
               fill = Shore),
           stat = "Identity",
           position = "dodge") +
  ylab("Precipitation (mm)") +
  xlab("Date") +
  scale_x_date(breaks = seq.Date(from = as.Date("2018-01-01"),
                                 to = as.Date("2019-07-01"),
                                 length.out = 4),
               labels = c("Jan\n2018","July\n2018","Jan\n2019","July\n2019")) +
  scale_fill_manual(values = carto_pal(12,"Safe")[c(2,9,10)]) +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.867,0.78),
        plot.margin = margin(0.1,0.15,0,0.1,"in"))
figure2

#### Save results ----
ggsave("plots/figure2.jpg",figure2,
       width = 4, height = 3, units = "in", dpi = 600)

#### Create recent rain sum
rollingPrecip = precip %>%
  mutate(Shore = case_when(Shore == "North" ~ "n",
                           Shore == "East" ~ "e",
                           Shore == "West" ~ "w")) %>%
  group_by(Date, Shore) %>%
  summarize(Precip_mm = mean(Precip_mm, na.rm = T)) %>%
  group_by(Shore) %>%
  mutate(rollPrecip = zoo::rollapply(Precip_mm, 
                                     width = 3, 
                                     FUN = "sum",
                                     align = "right",
                                     fill = NA))
print(rollingPrecip,n = 50)
rollingPrecip %>%
  ggplot(aes(x = Date, y = rollPrecip, group = Shore, col = Shore)) +
  geom_line()
write_csv(rollingPrecip,
          "data/rolling_precip.csv")

