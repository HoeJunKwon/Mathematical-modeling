library(readr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(corrplot)


## Data Load
male <- read_csv("Desktop/Male.csv")
female <- read_csv("Desktop/Female.csv")
male_route_infection <- read_csv("Desktop/Male_route_infection.csv")
female_route_infection <- read_csv("Desktop/Female_route_infection.csv")


## Data Summary

glimpse(male)
summary(male)

glimpse(female)
summary(female)

glimpse(male_route_infection)
summary(male_route_infection)

glimpse(female_route_infection)
summary(female_route_infection)


## Data Preprocessing

new_data1<- as.data.frame(apply(male[,-1],2,sum))
new_data2<- as.data.frame(apply(female[,-1],2,sum))


new_data<-cbind(new_data1[,1],new_data2[,1])
colnames(new_data) <- c("Male","Female")
Year <- seq(1985,2018,1)
new_data <- as.data.frame(cbind(Year,new_data))


cumMale<-cumsum(new_data$Male)
cumFemale<-cumsum(new_data$Female)
new_data_sum <-  as.data.frame(cbind(Year,cumMale,cumFemale))


vec<- c(apply(new_data_sum[,-1],1,sum))
vdc<-as.data.frame(vec)
vdc<-vdc[c(16:25),]


melt_data <- melt(data=new_data,
                  id.vars =c("Year"),
                  measure.vars = c("Male","Female")
)

melt_data_sum <- melt(data=new_data_sum,
                  id.vars =c("Year"),
                  measure.vars = c("cumMale","cumFemale")
)



melt_male_data <- melt(data = male_route_infection,
                       id.vars =c("Year"),
                       measure.vars = c("Heterosexual","Homosexual","Vertical_Transmission","Drug_Inject","Blood_Transfusion","No_Respose")
)

melt_female_data <- melt(data = female_route_infection,
                         id.vars =c("Year"),
                         measure.vars = c("Heterosexual","Homosexual","Vertical_Transmission","Drug_Inject","Blood_Transfusion","No_Respose")
)

# names(melt_male_data) <- c("Year","Route_of_Infection","Population")
# names(melt_female_data) <- c("Year","Route_of_Infection","Population")




# Data Visualization

ggplot(melt_data,aes(x=Year, y=value, group_by(variable),color=variable)) + geom_line() + geom_point()

ggplot(melt_data_sum,aes(x=Year, y=value, group_by(variable),color=variable)) + geom_line() + geom_point()

ggplot(melt_male_data,aes(x=Year, y=value, group_by(variable),color=variable)) + geom_line() + geom_point()

ggplot(melt_female_data,aes(x=Year, y=value, group_by(variable),color=variable)) + geom_line() + geom_point()

## Route of infection Coorelation plot

Cor_matrix = cor(male_route_infection[,2:7])
corrplot(Cor_matrix , method = "color",addgrid.col = "darkgray",addCoef.col = "white",
         tl.col = "indianred4",rect.col = "black",col = colorRampPalette(c("darkred","white","midnightblue"))(100))

## annual Route of infection box plot by sex

ggplot(melt_male_data,aes(x=variable, y=value,color=variable)) + geom_boxplot()
ggplot(melt_female_data,aes(x=variable, y=value,color=variable)) + geom_boxplot()













