load("C:/Document/Serieux/Travail/Package/comoR/inst/fit_plot_Neurips.RData")
df = data.frame(x=x,y=y, factor=as.factor(factor))


P1 <- ggplot(df, aes(x,y, col=factor))+geom_point(size=2)+
  theme_minimal()+theme( axis.text.y=element_blank(),

                         axis.ticks.y=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())


P1
df <- data.frame(x=x,y=y, L= fit_custom$L_pm[,1])
P11 <- ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("cEBMF")+theme_minimal()+theme( axis.text.y=element_blank(),

                                          axis.ticks.y=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank())


df <- data.frame(x=x,y=y, L= fit_default$L_pm[,1])
P12 <- ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("flashier")+theme_minimal()+theme( axis.text.y=element_blank(),

                                             axis.ticks.y=element_blank(),
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank())




df <- data.frame(x=x,y=y, L= fit_custom$L_pm[,2])

P21 <-ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("cEBMF")+theme_minimal()+theme( axis.text.y=element_blank(),

                                          axis.ticks.y=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank())


df <- data.frame(x=x,y=y, L= fit_default$L_pm[,2])
P22 <-ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("flashier")+theme_minimal()+
  theme( axis.text.y=element_blank(),

         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())


df <- data.frame(x=x,y=y, L= fit_custom$L_pm[,3])
P31 <-ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("cEBMF")+theme_minimal()+theme( axis.text.y=element_blank(),

                                          axis.ticks.y=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank())

df <- data.frame(x=x,y=y, L= fit_default$L_pm[,3])
P32 <-ggplot(df, aes ( x,y, col = abs(L)))+
  geom_point(size=2)+
  scale_color_continuous(type = "viridis")+
  geom_hline(yintercept = 0.33)+
  geom_hline(yintercept = 0.66)+
  geom_vline(xintercept = 0.66)+
  geom_vline(xintercept = 0.33)+
  ggtitle("flashier")+theme_minimal()+theme( axis.text.y=element_blank(),

                                             axis.ticks.y=element_blank(),
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank())

P1

gridExtra::grid.arrange(P11, P12, P31, P32, P21, P22, ncol=2)
