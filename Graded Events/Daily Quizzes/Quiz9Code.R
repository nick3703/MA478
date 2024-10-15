feedback <- read.csv("Course End Feedback Data.csv")

feedback_df <- data.frame(Section=ifelse(feedback$section=="MA206    H2 1 - CLARK, N",
                                         "HHour","IHour"),
                          HelpfulPeers= ifelse(feedback$My.fellow.students.contributed.to.my.learning.in.this.course. %in%
                                                 c("Strongly Agree","Agree"),"Yes","No"),
                          DidCourseHelpCreativity = feedback$n.this.course..my.ability.to.think.creatively.and.take.intellectual.risks.increased.)



fit.mod <- vglm(DidCourseHelpCreativity~Section,data=feedback_df,
                family=cumulative(parallel=TRUE))


fit.mod2 <- vglm(DidCourseHelpCreativity~1,data=feedback_df,
                family=cumulative(parallel=TRUE))