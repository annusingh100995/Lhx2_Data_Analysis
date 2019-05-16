STAGES = c("2_4_cell_1","1Kcell","dome","shield","bud","28hpf","2dpf","5dpf")

count_yap1 = c(132.8961962726,216.5132400557,257.9468053558,168.6918314337,269.5576905604,189.9377180122,150.2697454328,111.9334816179)
sd_yap1 = c(1.245942693,4.2740169122,3.1785777995,4.5401962406,1.0651122268,1.790694581,3.3142369716,1.3502556838)
Yap1 = data.frame(STAGES,count_yap1,sd_yap1)
Yap1$v1 <- factor(Yap1$STAGES,levels = Yap1$STAGES)
p1 = ggplot(Yap1, aes(x = v1, y = count_yap1))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = count_yap1-sd_yap1, ymax = (count_yap1)+(sd_yap1),width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Yap1",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission")
p3

count_sp1 = c(175.4629237295,181.0004680002,193.5971551014,119.9190375359,29.4511028832,33.33313479,22.9698447892,11.6675943714)
sd_sp1 = c(1.1039740297	,1.4707281667	,1.073763166,	14.5017877862,2.0330175315,0.8479289568,0.2658703717,1.5288555495)
sp1 = data.frame(STAGES,count_sp1,sd_sp1)
sp1$v1 <- factor(sp1$STAGES,levels = sp1$STAGES)
p1 = ggplot(sp1, aes(x = v1, y = count_sp1))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = count_sp1 - sd_sp1, ymax = count_sp1 + sd_sp1,width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Sp1",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission")
p3

count_ctnn1 = c(531.0571178543,446.4950365579,871.8301905091,932.7494957831,1423.3204076694,653.0773531979,463.2307721174,258.7936661888)
sd_ctnn1 = c(2.8488625941,11.24860831,14.430731334,	56.8066839857,15.5502743819,9.6217218357,8.3160633938,11.5630602362)
ctnn1 = data.frame(STAGES,count_ctnn1,sd_ctnn1)
ctnn1$v1 <- factor(ctnn1$STAGES,levels = ctnn1$STAGES)
p1 = ggplot(ctnn1, aes(x = v1, y = count_ctnn1))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = count_ctnn1 - sd_ctnn1, ymax = count_ctnn1 + sd_ctnn1,width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Ctnn1",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission")
p3

count_taz = c(0.3819432449,1.4462656793,3.6544709857,5.2567928118,11.5977975786,24.5115130409,17.8349190186,23.9663400097)
sd_taz = c(0.0949359338,0.1056778441,0.2775232571,0.2603203433,0.8863787902,0.2479902953,0.163632184,2.5040046614)
taz = data.frame(STAGES,count_taz,sd_taz)
taz$v1 <- factor(taz$STAGES,levels = taz$STAGES)
p1 = ggplot(taz, aes(x = v1, y = count_taz))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = count_taz - sd_taz, ymax = count_taz + sd_taz , width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Taz",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission")
p3


STAGES = c("2_4_cell_1","1Kcell","dome","shield","bud","28hpf","2dpf","5dpf")

count_yap1 = c(132.8961962726,216.5132400557,257.9468053558,168.6918314337,269.5576905604,189.9377180122,150.2697454328,111.9334816179)
sd_yap1 = c(1.245942693,4.2740169122,3.1785777995,4.5401962406,1.0651122268,1.790694581,3.3142369716,1.3502556838)
Yap1 = data.frame(STAGES,count_yap1,sd_yap1)
Yap1$v1 <- factor(Yap1$STAGES,levels = Yap1$STAGES)
p1 = ggplot(Yap1, aes(x = v1, y = count_yap1))+geom_bar(stat="identity", colour = "blue", fill = "blue")+geom_errorbar(aes(ymin = count_yap1-sd_yap1, ymax = (count_yap1)+(sd_yap1),width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Yap1",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission",text = element_text(size=40))
p3

count_sp1 = c(175.4629237295,181.0004680002,193.5971551014,119.9190375359,29.4511028832,33.33313479,22.9698447892,11.6675943714)
sd_sp1 = c(1.1039740297	,1.4707281667	,1.073763166,	14.5017877862,2.0330175315,0.8479289568,0.2658703717,1.5288555495)
sp1 = data.frame(STAGES,count_sp1,sd_sp1)
sp1$v1 <- factor(sp1$STAGES,levels = sp1$STAGES)
p1 = ggplot(sp1, aes(x = v1, y = count_sp1))+geom_bar(stat="identity", colour = "blue", fill = "blue")+geom_errorbar(aes(ymin = count_sp1 - sd_sp1, ymax = count_sp1 + sd_sp1,width = 0.2))+geom_line()
p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))
p3 = p2+ labs( title = "Sp1",x = "STAGES",y ="CPM",colour ="Cylinders",shape ="Transmission")
p3