      inv_chi = one/chi
      t1 = gtd(1,1)
      t2 = gtd(2,2)
      t4 = gtd(3,3)
      t6 = gtd(2,3)
      t7 = t6**2
      t9 = gtd(1,2)
      t10 = t9**2
      t12 = gtd(1,3)
      t16 = t12**2
      detgtd = t1*t2*t4-t1*t7-t10*t4+2*t9*t12*t6-t16*t2
      idetgtd = one/detgtd
      t2 = gtd(2,3)**2
      t14 = gtd(1,3)**2
      t22 = gtd(1,2)**2
      gtu(1,1) = idetgtd*(gtd(2,2)*gtd(3,3)-t2)
      gtu(1,2) = idetgtd*(-gtd(1,2)*gtd(3,3)+gtd(1,3)*gtd(2,3))
      gtu(1,3) = idetgtd*(gtd(1,2)*gtd(2,3)-gtd(1,3)*gtd(2,2))
      gtu(2,1) = gtu(1,2)
      gtu(2,2) = idetgtd*(gtd(1,1)*gtd(3,3)-t14)
      gtu(2,3) = idetgtd*(-gtd(1,1)*gtd(2,3)+gtd(1,2)*gtd(1,3))
      gtu(3,1) = gtu(1,3)
      gtu(3,2) = gtu(2,3)
      gtu(3,3) = idetgtd*(gtd(1,1)*gtd(2,2)-t22)
      t2 = gtu(1,2)*Atd(1,2)
      t3 = gtu(1,3)*Atd(1,3)
      t18 = gtu(2,3)*Atd(2,3)
      Atud(1,1) = gtu(1,1)*Atd(1,1)+t2+t3
      Atud(1,2) = gtu(1,1)*Atd(1,2)+gtu(1,2)*Atd(2,2)+gtu(1,3)*Atd(2,3)
      Atud(1,3) = gtu(1,1)*Atd(1,3)+gtu(1,2)*Atd(2,3)+gtu(1,3)*Atd(3,3)
      Atud(2,1) = gtu(1,2)*Atd(1,1)+gtu(2,2)*Atd(1,2)+gtu(2,3)*Atd(1,3)
      Atud(2,2) = t2+gtu(2,2)*Atd(2,2)+t18
      Atud(2,3) = gtu(1,2)*Atd(1,3)+gtu(2,2)*Atd(2,3)+gtu(2,3)*Atd(3,3)
      Atud(3,1) = gtu(1,3)*Atd(1,1)+gtu(2,3)*Atd(1,2)+gtu(3,3)*Atd(1,3)
      Atud(3,2) = gtu(1,3)*Atd(1,2)+gtu(2,3)*Atd(2,2)+gtu(3,3)*Atd(2,3)
      Atud(3,3) = t3+t18+gtu(3,3)*Atd(3,3)
      Atu(1,1) = Atud(1,1)*gtu(1,1)+Atud(1,2)*gtu(1,2)+Atud(1,3)*gtu(1,3
     #)
      Atu(1,2) = Atud(1,1)*gtu(1,2)+Atud(1,2)*gtu(2,2)+Atud(1,3)*gtu(2,3
     #)
      Atu(1,3) = Atud(1,1)*gtu(1,3)+Atud(1,2)*gtu(2,3)+Atud(1,3)*gtu(3,3
     #)
      Atu(2,1) = Atu(1,2)
      Atu(2,2) = Atud(2,1)*gtu(1,2)+Atud(2,2)*gtu(2,2)+Atud(2,3)*gtu(2,3
     #)
      Atu(2,3) = Atud(2,1)*gtu(1,3)+Atud(2,2)*gtu(2,3)+Atud(2,3)*gtu(3,3
     #)
      Atu(3,1) = Atu(1,3)
      Atu(3,2) = Atu(2,3)
      Atu(3,3) = Atud(3,1)*gtu(1,3)+Atud(3,2)*gtu(2,3)+Atud(3,3)*gtu(3,3
     #)
      Ctd(1,1,1) = d_gtd(1,1,1)
      Ctd(1,1,2) = d_gtd(2,1,1)
      Ctd(1,1,3) = d_gtd(3,1,1)
      Ctd(1,2,1) = Ctd(1,1,2)
      Ctd(1,2,2) = 2*d_gtd(2,1,2)-d_gtd(1,2,2)
      Ctd(1,2,3) = d_gtd(2,1,3)+d_gtd(3,1,2)-d_gtd(1,2,3)
      Ctd(1,3,1) = Ctd(1,1,3)
      Ctd(1,3,2) = Ctd(1,2,3)
      Ctd(1,3,3) = 2*d_gtd(3,1,3)-d_gtd(1,3,3)
      Ctd(2,1,1) = 2*d_gtd(1,1,2)-d_gtd(2,1,1)
      Ctd(2,1,2) = d_gtd(1,2,2)
      Ctd(2,1,3) = d_gtd(1,2,3)+d_gtd(3,1,2)-d_gtd(2,1,3)
      Ctd(2,2,1) = Ctd(2,1,2)
      Ctd(2,2,2) = d_gtd(2,2,2)
      Ctd(2,2,3) = d_gtd(3,2,2)
      Ctd(2,3,1) = Ctd(2,1,3)
      Ctd(2,3,2) = Ctd(2,2,3)
      Ctd(2,3,3) = 2*d_gtd(3,2,3)-d_gtd(2,3,3)
      Ctd(3,1,1) = 2*d_gtd(1,1,3)-d_gtd(3,1,1)
      Ctd(3,1,2) = d_gtd(1,2,3)+d_gtd(2,1,3)-d_gtd(3,1,2)
      Ctd(3,1,3) = d_gtd(1,3,3)
      Ctd(3,2,1) = Ctd(3,1,2)
      Ctd(3,2,2) = 2*d_gtd(2,2,3)-d_gtd(3,2,2)
      Ctd(3,2,3) = d_gtd(2,3,3)
      Ctd(3,3,1) = Ctd(3,1,3)
      Ctd(3,3,2) = Ctd(3,2,3)
      Ctd(3,3,3) = d_gtd(3,3,3)
      Ct(1,1,1) = gtu(1,1)*Ctd(1,1,1)+gtu(1,2)*Ctd(2,1,1)+gtu(1,3)*Ctd(3
     #,1,1)
      Ct(1,1,2) = gtu(1,1)*Ctd(1,1,2)+gtu(1,2)*Ctd(2,1,2)+gtu(1,3)*Ctd(3
     #,1,2)
      Ct(1,1,3) = gtu(1,1)*Ctd(1,1,3)+gtu(1,2)*Ctd(2,1,3)+gtu(1,3)*Ctd(3
     #,1,3)
      Ct(1,2,1) = Ct(1,1,2)
      Ct(1,2,2) = gtu(1,1)*Ctd(1,2,2)+gtu(1,2)*Ctd(2,2,2)+gtu(1,3)*Ctd(3
     #,2,2)
      Ct(1,2,3) = gtu(1,1)*Ctd(1,2,3)+gtu(1,2)*Ctd(2,2,3)+gtu(1,3)*Ctd(3
     #,2,3)
      Ct(1,3,1) = Ct(1,1,3)
      Ct(1,3,2) = Ct(1,2,3)
      Ct(1,3,3) = gtu(1,1)*Ctd(1,3,3)+gtu(1,2)*Ctd(2,3,3)+gtu(1,3)*Ctd(3
     #,3,3)
      Ct(2,1,1) = gtu(1,2)*Ctd(1,1,1)+gtu(2,2)*Ctd(2,1,1)+gtu(2,3)*Ctd(3
     #,1,1)
      Ct(2,1,2) = gtu(1,2)*Ctd(1,1,2)+gtu(2,2)*Ctd(2,1,2)+gtu(2,3)*Ctd(3
     #,1,2)
      Ct(2,1,3) = gtu(1,2)*Ctd(1,1,3)+gtu(2,2)*Ctd(2,1,3)+gtu(2,3)*Ctd(3
     #,1,3)
      Ct(2,2,1) = Ct(2,1,2)
      Ct(2,2,2) = gtu(1,2)*Ctd(1,2,2)+gtu(2,2)*Ctd(2,2,2)+gtu(2,3)*Ctd(3
     #,2,2)
      Ct(2,2,3) = gtu(1,2)*Ctd(1,2,3)+gtu(2,2)*Ctd(2,2,3)+gtu(2,3)*Ctd(3
     #,2,3)
      Ct(2,3,1) = Ct(2,1,3)
      Ct(2,3,2) = Ct(2,2,3)
      Ct(2,3,3) = gtu(1,2)*Ctd(1,3,3)+gtu(2,2)*Ctd(2,3,3)+gtu(2,3)*Ctd(3
     #,3,3)
      Ct(3,1,1) = gtu(1,3)*Ctd(1,1,1)+gtu(2,3)*Ctd(2,1,1)+gtu(3,3)*Ctd(3
     #,1,1)
      Ct(3,1,2) = gtu(1,3)*Ctd(1,1,2)+gtu(2,3)*Ctd(2,1,2)+gtu(3,3)*Ctd(3
     #,1,2)
      Ct(3,1,3) = gtu(1,3)*Ctd(1,1,3)+gtu(2,3)*Ctd(2,1,3)+gtu(3,3)*Ctd(3
     #,1,3)
      Ct(3,2,1) = Ct(3,1,2)
      Ct(3,2,2) = gtu(1,3)*Ctd(1,2,2)+gtu(2,3)*Ctd(2,2,2)+gtu(3,3)*Ctd(3
     #,2,2)
      Ct(3,2,3) = gtu(1,3)*Ctd(1,2,3)+gtu(2,3)*Ctd(2,2,3)+gtu(3,3)*Ctd(3
     #,2,3)
      Ct(3,3,1) = Ct(3,1,3)
      Ct(3,3,2) = Ct(3,2,3)
      Ct(3,3,3) = gtu(1,3)*Ctd(1,3,3)+gtu(2,3)*Ctd(2,3,3)+gtu(3,3)*Ctd(3
     #,3,3)
      div_Beta = d_Betau(1,1)+d_Betau(2,2)+d_Betau(3,3)
      d_div_Beta(1) = dd_Betau(1,1,1)+dd_Betau(1,2,2)+dd_Betau(1,3,3)
      d_div_Beta(2) = dd_Betau(1,2,1)+dd_Betau(2,2,2)+dd_Betau(2,3,3)
      d_div_Beta(3) = dd_Betau(1,3,1)+dd_Betau(2,3,2)+dd_Betau(3,3,3)
      CalGamt(1) = half*(gtu(1,1)*Ct(1,1,1)+gtu(2,2)*Ct(1,2,2)+gtu(3,3)*
     #Ct(1,3,3))+gtu(1,2)*Ct(1,1,2)+gtu(1,3)*Ct(1,1,3)+gtu(2,3)*Ct(1,2,3
     #)
      CalGamt(2) = half*(gtu(1,1)*Ct(2,1,1)+gtu(2,2)*Ct(2,2,2)+gtu(3,3)*
     #Ct(2,3,3))+gtu(1,2)*Ct(2,1,2)+gtu(1,3)*Ct(2,1,3)+gtu(2,3)*Ct(2,2,3
     #)
      CalGamt(3) = half*(gtu(1,1)*Ct(3,1,1)+gtu(2,2)*Ct(3,2,2)+gtu(3,3)*
     #Ct(3,3,3))+gtu(1,2)*Ct(3,1,2)+gtu(1,3)*Ct(3,1,3)+gtu(2,3)*Ct(3,2,3
     #)
      t2 = d_chi(1)**2
      t29 = d_chi(2)**2
      t47 = d_chi(3)**2
      Rpd_1(1,1) = half*dd_chi(1,1)-fourth*(t2*inv_chi+Ct(1,1,1)*d_chi(1
     #)+Ct(2,1,1)*d_chi(2)+Ct(3,1,1)*d_chi(3))
      Rpd_1(1,2) = half*dd_chi(1,2)-fourth*(d_chi(1)*d_chi(2)*inv_chi+Ct
     #(1,1,2)*d_chi(1)+Ct(2,1,2)*d_chi(2)+Ct(3,1,2)*d_chi(3))
      Rpd_1(1,3) = half*dd_chi(1,3)-fourth*(d_chi(1)*d_chi(3)*inv_chi+Ct
     #(1,1,3)*d_chi(1)+Ct(2,1,3)*d_chi(2)+Ct(3,1,3)*d_chi(3))
      Rpd_1(2,1) = Rpd_1(1,2)
      Rpd_1(2,2) = half*dd_chi(2,2)-fourth*(t29*inv_chi+Ct(1,2,2)*d_chi(
     #1)+Ct(2,2,2)*d_chi(2)+Ct(3,2,2)*d_chi(3))
      Rpd_1(2,3) = half*dd_chi(2,3)-fourth*(d_chi(2)*d_chi(3)*inv_chi+Ct
     #(1,2,3)*d_chi(1)+Ct(2,2,3)*d_chi(2)+Ct(3,2,3)*d_chi(3))
      Rpd_1(3,1) = Rpd_1(1,3)
      Rpd_1(3,2) = Rpd_1(2,3)
      Rpd_1(3,3) = half*dd_chi(3,3)-fourth*(t47*inv_chi+Ct(1,3,3)*d_chi(
     #1)+Ct(2,3,3)*d_chi(2)+Ct(3,3,3)*d_chi(3))

!---------------------------------------------------------------------
      t2 = d_chi(1)**2
      t17 = d_chi(2)**2
      t22 = d_chi(3)**2
      t28 = d_chi(1)*inv_chi
      t32 = threehalves*d_chi(3)
      t42 = CalGamt(1)*d_chi(1)+CalGamt(2)*d_chi(2)+CalGamt(3)*d_chi(3)+
     #gtu(1,1)*(threehalves*t2*inv_chi-dd_chi(1,1))+gtu(2,2)*(threehalve
     #s*t17*inv_chi-dd_chi(2,2))+gtu(3,3)*(threehalves*t22*inv_chi-dd_ch
     #i(3,3))+two*(gtu(1,2)*(threehalves*d_chi(2)*t28-dd_chi(1,2))+gtu(1
     #,3)*(t32*t28-dd_chi(1,3))+gtu(2,3)*(t32*d_chi(2)*inv_chi-dd_chi(2,
     #3)))
      Rpd(1,1) = half*dd_chi(1,1)-fourth*(t2*inv_chi+Ct(1,1,1)*d_chi(1)+
     #Ct(2,1,1)*d_chi(2)+Ct(3,1,1)*d_chi(3))-half*gtd(1,1)*t42
      Rpd(1,2) = half*dd_chi(1,2)-fourth*(d_chi(1)*d_chi(2)*inv_chi+Ct(1
     #,1,2)*d_chi(1)+Ct(2,1,2)*d_chi(2)+Ct(3,1,2)*d_chi(3))-half*gtd(1,2
     #)*t42
      Rpd(1,3) = half*dd_chi(1,3)-fourth*(d_chi(1)*d_chi(3)*inv_chi+Ct(1
     #,1,3)*d_chi(1)+Ct(2,1,3)*d_chi(2)+Ct(3,1,3)*d_chi(3))-half*gtd(1,3
     #)*t42
      Rpd(2,1) = Rpd(1,2)
      Rpd(2,2) = half*dd_chi(2,2)-fourth*(t17*inv_chi+Ct(1,2,2)*d_chi(1)
     #+Ct(2,2,2)*d_chi(2)+Ct(3,2,2)*d_chi(3))-half*gtd(2,2)*t42
      Rpd(2,3) = half*dd_chi(2,3)-fourth*(d_chi(2)*d_chi(3)*inv_chi+Ct(1
     #,2,3)*d_chi(1)+Ct(2,2,3)*d_chi(2)+Ct(3,2,3)*d_chi(3))-half*gtd(2,3
     #)*t42
      Rpd(3,1) = Rpd(1,3)
      Rpd(3,2) = Rpd(2,3)
      Rpd(3,3) = half*dd_chi(3,3)-fourth*(t22*inv_chi+Ct(1,3,3)*d_chi(1)
     #+Ct(2,3,3)*d_chi(2)+Ct(3,3,3)*d_chi(3))-half*gtd(3,3)*t42


!---------------------------------------------------------------------




      t17 = Ct(1,1,2)*Ctd(1,1,1)
      t19 = Ct(1,1,1)*Ctd(1,1,2)
      t20 = Ct(2,1,2)*Ctd(1,1,2)
      t22 = Ct(2,1,1)*Ctd(2,1,2)
      t23 = Ct(3,1,2)*Ctd(1,1,3)
      t25 = Ct(3,1,1)*Ctd(3,1,2)
      t28 = Ct(1,1,3)*Ctd(1,1,1)
      t30 = Ct(1,1,1)*Ctd(1,1,3)
      t31 = Ct(2,1,3)*Ctd(1,1,2)
      t33 = Ct(2,1,1)*Ctd(2,1,3)
      t34 = Ct(3,1,3)*Ctd(1,1,3)
      t36 = Ct(3,1,1)*Ctd(3,1,3)
      t48 = Ct(1,1,2)*Ctd(1,1,2)
      t50 = Ct(2,1,2)*Ctd(1,2,2)
      t52 = Ct(2,1,2)*Ctd(2,1,2)
      t53 = Ct(3,1,2)*Ctd(1,2,3)
      t55 = Ct(3,1,2)*Ctd(3,1,2)
      t58 = Ct(1,1,3)*Ctd(1,1,2)
      t60 = Ct(1,1,2)*Ctd(1,1,3)
      t61 = Ct(2,1,3)*Ctd(1,2,2)
      t63 = Ct(2,1,2)*Ctd(2,1,3)
      t64 = Ct(3,1,3)*Ctd(1,2,3)
      t66 = Ct(3,1,2)*Ctd(3,1,3)
      t79 = Ct(2,1,2)*Ctd(1,2,3)
      t81 = Ct(2,1,3)*Ctd(2,1,2)
      t82 = Ct(3,1,2)*Ctd(1,3,3)
      t84 = Ct(3,1,3)*Ctd(3,1,2)
      t87 = Ct(1,1,3)*Ctd(1,1,3)
      t89 = Ct(2,1,3)*Ctd(1,2,3)
      t91 = Ct(2,1,3)*Ctd(2,1,3)
      t92 = Ct(3,1,3)*Ctd(1,3,3)
      t94 = Ct(3,1,3)*Ctd(3,1,3)
      t97 = 2*CalGamt(1)*Ctd(1,1,1)+2*CalGamt(2)*Ctd(1,1,2)+2*CalGamt(3)
     #*Ctd(1,1,3)+gtu(1,1)*(3*Ct(1,1,1)*Ctd(1,1,1)+2*Ct(2,1,1)*Ctd(1,1,2
     #)+Ct(2,1,1)*Ctd(2,1,1)+2*Ct(3,1,1)*Ctd(1,1,3)+Ct(3,1,1)*Ctd(3,1,1)
     #)+gtu(1,2)*(2*t17+t19+2*t20+t22+2*t23+t25)+gtu(1,3)*(2*t28+t30+2*t
     #31+t33+2*t34+t36)+gtu(1,2)*(2*t19+t17+2*Ct(2,1,1)*Ctd(1,2,2)+Ct(2,
     #1,2)*Ctd(2,1,1)+2*Ct(3,1,1)*Ctd(1,2,3)+Ct(3,1,2)*Ctd(3,1,1))+gtu(2
     #,2)*(3*t48+2*t50+t52+2*t53+t55)+gtu(2,3)*(2*t58+t60+2*t61+t63+2*t6
     #4+t66)+gtu(1,3)*(2*t30+t28+2*Ct(2,1,1)*Ctd(1,2,3)+Ct(2,1,3)*Ctd(2,
     #1,1)+2*Ct(3,1,1)*Ctd(1,3,3)+Ct(3,1,3)*Ctd(3,1,1))+gtu(2,3)*(2*t60+
     #t58+2*t79+t81+2*t82+t84)+gtu(3,3)*(3*t87+2*t89+t91+2*t92+t94)
      t125 = Ct(1,1,2)*Ctd(2,1,1)
      t129 = Ct(2,1,1)*Ctd(2,2,2)
      t130 = Ct(3,1,2)*Ctd(2,1,3)
      t135 = Ct(1,1,3)*Ctd(2,1,1)
      t136 = Ct(1,2,3)*Ctd(1,1,1)
      t137 = Ct(1,1,1)*Ctd(1,2,3)
      t138 = Ct(2,2,3)*Ctd(1,1,2)
      t139 = Ct(2,1,1)*Ctd(2,2,3)
      t140 = Ct(3,1,3)*Ctd(2,1,3)
      t141 = Ct(3,2,3)*Ctd(1,1,3)
      t142 = Ct(3,1,1)*Ctd(3,2,3)
      t150 = Ct(1,1,2)*Ctd(2,1,2)
      t151 = Ct(1,2,2)*Ctd(1,1,2)
      t152 = Ct(1,1,2)*Ctd(1,2,2)
      t153 = Ct(2,1,2)*Ctd(2,2,2)
      t154 = 2*t153
      t156 = Ct(3,1,2)*Ctd(2,2,3)
      t158 = Ct(3,1,2)*Ctd(3,2,2)
      t161 = Ct(1,1,3)*Ctd(2,1,2)
      t162 = Ct(1,2,3)*Ctd(1,1,2)
      t163 = Ct(1,1,2)*Ctd(1,2,3)
      t164 = Ct(2,1,3)*Ctd(2,2,2)
      t165 = Ct(2,2,3)*Ctd(1,2,2)
      t166 = Ct(2,1,2)*Ctd(2,2,3)
      t167 = Ct(3,1,3)*Ctd(2,2,3)
      t168 = Ct(3,2,3)*Ctd(1,2,3)
      t169 = Ct(3,1,2)*Ctd(3,2,3)
      t176 = Ct(1,1,2)*Ctd(2,1,3)
      t177 = Ct(1,2,2)*Ctd(1,1,3)
      t180 = Ct(3,1,2)*Ctd(2,3,3)
      t185 = Ct(1,1,3)*Ctd(2,1,3)
      t186 = Ct(1,2,3)*Ctd(1,1,3)
      t187 = Ct(1,1,3)*Ctd(1,2,3)
      t188 = Ct(2,1,3)*Ctd(2,2,3)
      t190 = Ct(2,2,3)*Ctd(1,2,3)
      t191 = Ct(3,1,3)*Ctd(2,3,3)
      t192 = Ct(3,2,3)*Ctd(1,3,3)
      t193 = Ct(3,1,3)*Ctd(3,2,3)
      s1 = CalGamt(1)*(Ctd(1,1,2)+Ctd(2,1,1))+CalGamt(2)*(Ctd(1,2,2)+Ctd
     #(2,1,2))+CalGamt(3)*(Ctd(1,2,3)+Ctd(2,1,3))+gtu(1,1)*(Ct(1,1,1)*Ct
     #d(2,1,1)+t17+t19+2*t22+t20+Ct(3,1,1)*Ctd(2,1,3)+t23+t25)+gtu(1,2)*
     #(t125+Ct(1,2,2)*Ctd(1,1,1)+Ct(1,1,1)*Ctd(1,2,2)+t52+Ct(2,2,2)*Ctd(
     #1,1,2)+t129+t130+Ct(3,2,2)*Ctd(1,1,3)+Ct(3,1,1)*Ctd(3,2,2))+gtu(1,
     #3)*(t135+t136+t137+t81+t138+t139+t140+t141+t142)
      t196 = s1+gtu(1,2)*(Ct(1,1,1)*Ctd(2,1,2)+2*t48+t129+t50+t52+Ct(3,1
     #,1)*Ctd(2,2,3)+t53+t55)+gtu(2,2)*(t150+t151+t152+t154+Ct(2,2,2)*Ct
     #d(1,2,2)+t156+Ct(3,2,2)*Ctd(1,2,3)+t158)+gtu(2,3)*(t161+t162+t163+
     #t164+t165+t166+t167+t168+t169)+gtu(1,3)*(Ct(1,1,1)*Ctd(2,1,3)+t60+
     #t58+t139+t79+t81+Ct(3,1,1)*Ctd(2,3,3)+t82+t84)+gtu(2,3)*(t176+t177
     #+Ct(1,1,3)*Ctd(1,2,2)+t166+Ct(2,2,2)*Ctd(1,2,3)+t164+t180+Ct(3,2,2
     #)*Ctd(1,3,3)+Ct(3,1,3)*Ctd(3,2,2))+gtu(3,3)*(t185+t186+t187+2*t188
     #+t190+t191+t192+t193)
      t224 = Ct(1,1,2)*Ctd(3,1,1)
      t225 = Ct(2,1,2)*Ctd(3,1,2)
      t228 = Ct(1,1,3)*Ctd(3,1,1)
      t231 = Ct(2,1,3)*Ctd(3,1,2)
      t235 = Ct(3,1,1)*Ctd(3,3,3)
      t242 = Ct(1,1,2)*Ctd(3,1,2)
      t243 = Ct(2,1,2)*Ctd(3,2,2)
      t247 = Ct(1,1,3)*Ctd(3,1,2)
      t249 = Ct(1,1,2)*Ctd(1,3,3)
      t250 = Ct(2,1,3)*Ctd(3,2,2)
      t252 = Ct(2,1,2)*Ctd(2,3,3)
      t254 = Ct(3,1,2)*Ctd(3,3,3)
      t262 = Ct(1,1,2)*Ctd(3,1,3)
      t263 = Ct(2,1,2)*Ctd(3,2,3)
      t266 = Ct(1,1,3)*Ctd(3,1,3)
      t267 = Ct(1,3,3)*Ctd(1,1,3)
      t268 = Ct(1,1,3)*Ctd(1,3,3)
      t269 = Ct(2,1,3)*Ctd(3,2,3)
      t271 = Ct(2,1,3)*Ctd(2,3,3)
      t272 = Ct(3,1,3)*Ctd(3,3,3)
      t273 = 2*t272
      s1 = CalGamt(1)*(Ctd(1,1,3)+Ctd(3,1,1))+CalGamt(2)*(Ctd(1,2,3)+Ctd
     #(3,1,2))+CalGamt(3)*(Ctd(1,3,3)+Ctd(3,1,3))+gtu(1,1)*(Ct(1,1,1)*Ct
     #d(3,1,1)+t28+t30+Ct(2,1,1)*Ctd(3,1,2)+t31+t33+2*t36+t34)+gtu(1,2)*
     #(t224+t136+t137+t225+t138+t139+t66+t141+t142)+gtu(1,3)*(t228+Ct(1,
     #3,3)*Ctd(1,1,1)+Ct(1,1,1)*Ctd(1,3,3)+t231+Ct(2,3,3)*Ctd(1,1,2)+Ct(
     #2,1,1)*Ctd(2,3,3)+t94+Ct(3,3,3)*Ctd(1,1,3)+t235)
      t277 = s1+gtu(1,2)*(Ct(1,1,1)*Ctd(3,1,2)+t58+t60+Ct(2,1,1)*Ctd(3,2
     #,2)+t61+t63+t142+t64+t66)+gtu(2,2)*(t242+t162+t163+t243+t165+t166+
     #2*t169+t168)+gtu(2,3)*(t247+Ct(1,3,3)*Ctd(1,1,2)+t249+t250+Ct(2,3,
     #3)*Ctd(1,2,2)+t252+t193+Ct(3,3,3)*Ctd(1,2,3)+t254)+gtu(1,3)*(Ct(1,
     #1,1)*Ctd(3,1,3)+2*t87+Ct(2,1,1)*Ctd(3,2,3)+t89+t91+t235+t92+t94)+g
     #tu(2,3)*(t262+t186+t187+t263+t190+t188+t254+t192+t193)+gtu(3,3)*(t
     #266+t267+t268+t269+Ct(2,3,3)*Ctd(1,2,3)+t271+t273+Ct(3,3,3)*Ctd(1,
     #3,3))
      t307 = Ct(2,2,2)*Ctd(2,1,2)
      t313 = Ct(1,2,3)*Ctd(2,1,1)
      t315 = Ct(2,2,3)*Ctd(2,1,2)
      t317 = Ct(3,2,3)*Ctd(2,1,3)
      t336 = Ct(1,2,3)*Ctd(2,1,2)
      t338 = Ct(1,2,2)*Ctd(1,2,3)
      t339 = Ct(2,2,3)*Ctd(2,2,2)
      t341 = Ct(2,2,2)*Ctd(2,2,3)
      t342 = Ct(3,2,3)*Ctd(2,2,3)
      t344 = Ct(3,2,2)*Ctd(3,2,3)
      t362 = Ct(1,2,3)*Ctd(2,1,3)
      t364 = Ct(1,2,3)*Ctd(1,2,3)
      t365 = Ct(2,2,3)*Ctd(2,2,3)
      t367 = Ct(3,2,3)*Ctd(2,3,3)
      t369 = Ct(3,2,3)*Ctd(3,2,3)
      s1 = 2*CalGamt(1)*Ctd(2,1,2)+2*CalGamt(2)*Ctd(2,2,2)+2*CalGamt(3)*
     #Ctd(2,2,3)+gtu(1,1)*(2*t125+t48+3*t52+2*t130+t55)+gtu(1,2)*(2*Ct(1
     #,2,2)*Ctd(2,1,1)+t152+2*t307+t153+2*Ct(3,2,2)*Ctd(2,1,3)+t158)+gtu
     #(1,3)*(2*t313+t163+2*t315+t166+2*t317+t169)
      t372 = s1+gtu(1,2)*(2*t150+t151+t154+t307+2*t156+Ct(3,2,2)*Ctd(3,1
     #,2))+gtu(2,2)*(2*Ct(1,2,2)*Ctd(2,1,2)+Ct(1,2,2)*Ctd(1,2,2)+3*Ct(2,
     #2,2)*Ctd(2,2,2)+2*Ct(3,2,2)*Ctd(2,2,3)+Ct(3,2,2)*Ctd(3,2,2))+gtu(2
     #,3)*(2*t336+t338+2*t339+t341+2*t342+t344)+gtu(1,3)*(2*t176+t162+2*
     #t166+t315+2*t180+Ct(3,2,3)*Ctd(3,1,2))+gtu(2,3)*(2*Ct(1,2,2)*Ctd(2
     #,1,3)+Ct(1,2,3)*Ctd(1,2,2)+2*t341+t339+2*Ct(3,2,2)*Ctd(2,3,3)+Ct(3
     #,2,3)*Ctd(3,2,2))+gtu(3,3)*(2*t362+t364+3*t365+2*t367+t369)
      t400 = Ct(3,2,2)*Ctd(3,1,3)
      t403 = Ct(1,2,3)*Ctd(3,1,1)
      t405 = Ct(2,2,3)*Ctd(3,1,2)
      t407 = Ct(3,2,3)*Ctd(3,1,3)
      t419 = Ct(1,2,3)*Ctd(3,1,2)
      t422 = Ct(2,2,3)*Ctd(3,2,2)
      t426 = Ct(3,2,2)*Ctd(3,3,3)
      t429 = Ct(2,2,3)*Ctd(2,1,3)
      t437 = Ct(1,2,3)*Ctd(3,1,3)
      t439 = Ct(1,2,3)*Ctd(1,3,3)
      t440 = Ct(2,2,3)*Ctd(3,2,3)
      t441 = Ct(2,3,3)*Ctd(2,2,3)
      t442 = Ct(2,2,3)*Ctd(2,3,3)
      t443 = Ct(3,2,3)*Ctd(3,3,3)
      t444 = 2*t443
      s1 = CalGamt(1)*(Ctd(2,1,3)+Ctd(3,1,2))+CalGamt(2)*(Ctd(2,2,3)+Ctd
     #(3,2,2))+CalGamt(3)*(Ctd(2,3,3)+Ctd(3,2,3))+gtu(1,1)*(t224+t135+t6
     #0+t225+t81+t63+2*t66+t140)+gtu(1,2)*(Ct(1,2,2)*Ctd(3,1,1)+t313+t16
     #3+Ct(2,2,2)*Ctd(3,1,2)+t315+t166+t400+t317+t169)+gtu(1,3)*(t403+Ct
     #(1,3,3)*Ctd(2,1,1)+t249+t405+Ct(2,3,3)*Ctd(2,1,2)+t252+t407+Ct(3,3
     #,3)*Ctd(2,1,3)+t254)
      t448 = s1+gtu(1,2)*(t242+t161+t177+t243+t164+Ct(2,2,2)*Ctd(2,1,3)+
     #t169+t167+t400)+gtu(2,2)*(Ct(1,2,2)*Ctd(3,1,2)+t336+t338+Ct(2,2,2)
     #*Ctd(3,2,2)+t339+t341+2*t344+t342)+gtu(2,3)*(t419+Ct(1,3,3)*Ctd(2,
     #1,2)+Ct(1,2,2)*Ctd(1,3,3)+t422+Ct(2,3,3)*Ctd(2,2,2)+Ct(2,2,2)*Ctd(
     #2,3,3)+t369+Ct(3,3,3)*Ctd(2,2,3)+t426)+gtu(1,3)*(t262+t185+t186+t2
     #63+t188+t429+t254+t191+t407)+gtu(2,3)*(Ct(1,2,2)*Ctd(3,1,3)+t362+t
     #364+Ct(2,2,2)*Ctd(3,2,3)+2*t365+t426+t367+t369)+gtu(3,3)*(t437+Ct(
     #1,3,3)*Ctd(2,1,3)+t439+t440+t441+t442+t444+Ct(3,3,3)*Ctd(2,3,3))
      t485 = Ct(3,3,3)*Ctd(3,1,3)
      t503 = Ct(3,3,3)*Ctd(3,2,3)
      t527 = 2*CalGamt(1)*Ctd(3,1,3)+2*CalGamt(2)*Ctd(3,2,3)+2*CalGamt(3
     #)*Ctd(3,3,3)+gtu(1,1)*(2*t228+t87+2*t231+t91+3*t94)+gtu(1,2)*(2*t4
     #03+t187+2*t405+t188+2*t407+t193)+gtu(1,3)*(2*Ct(1,3,3)*Ctd(3,1,1)+
     #t268+2*Ct(2,3,3)*Ctd(3,1,2)+t271+2*t485+t272)+gtu(1,2)*(2*t247+t18
     #6+2*t250+t429+2*t193+t407)+gtu(2,2)*(2*t419+t364+2*t422+t365+3*t36
     #9)+gtu(2,3)*(2*Ct(1,3,3)*Ctd(3,1,2)+t439+2*Ct(2,3,3)*Ctd(3,2,2)+t4
     #42+2*t503+t443)+gtu(1,3)*(2*t266+t267+2*t269+Ct(2,3,3)*Ctd(2,1,3)+
     #t273+t485)+gtu(2,3)*(2*t437+Ct(1,3,3)*Ctd(1,2,3)+2*t440+t441+t444+
     #t503)+gtu(3,3)*(2*Ct(1,3,3)*Ctd(3,1,3)+Ct(1,3,3)*Ctd(1,3,3)+2*Ct(2
     #,3,3)*Ctd(3,2,3)+Ct(2,3,3)*Ctd(2,3,3)+3*Ct(3,3,3)*Ctd(3,3,3))
      Rtd(1,1) = fourth*t97-half*(gtu(1,1)*dd_gtd(1,1,1,1)+gtu(2,2)*dd_g
     #td(2,2,1,1)+gtu(3,3)*dd_gtd(3,3,1,1)-2*gtd(1,1)*d_Gamt(1,1)-2*gtd(
     #1,2)*d_Gamt(1,2)-2*gtd(1,3)*d_Gamt(1,3))-gtu(1,2)*dd_gtd(1,2,1,1)-
     #gtu(1,3)*dd_gtd(1,3,1,1)-gtu(2,3)*dd_gtd(2,3,1,1)
      Rtd(1,2) = fourth*t196-half*(gtu(1,1)*dd_gtd(1,1,1,2)+gtu(2,2)*dd_
     #gtd(2,2,1,2)+gtu(3,3)*dd_gtd(3,3,1,2)-gtd(1,1)*d_Gamt(2,1)-gtd(1,2
     #)*d_Gamt(1,1)-gtd(1,2)*d_Gamt(2,2)-gtd(2,2)*d_Gamt(1,2)-gtd(1,3)*d
     #_Gamt(2,3)-gtd(2,3)*d_Gamt(1,3))-gtu(1,2)*dd_gtd(1,2,1,2)-gtu(1,3)
     #*dd_gtd(1,3,1,2)-gtu(2,3)*dd_gtd(2,3,1,2)
      Rtd(1,3) = fourth*t277-half*(gtu(1,1)*dd_gtd(1,1,1,3)+gtu(2,2)*dd_
     #gtd(2,2,1,3)+gtu(3,3)*dd_gtd(3,3,1,3)-gtd(1,1)*d_Gamt(3,1)-gtd(1,3
     #)*d_Gamt(1,1)-gtd(1,2)*d_Gamt(3,2)-gtd(2,3)*d_Gamt(1,2)-gtd(1,3)*d
     #_Gamt(3,3)-gtd(3,3)*d_Gamt(1,3))-gtu(1,2)*dd_gtd(1,2,1,3)-gtu(1,3)
     #*dd_gtd(1,3,1,3)-gtu(2,3)*dd_gtd(2,3,1,3)
      Rtd(2,1) = Rtd(1,2)
      Rtd(2,2) = fourth*t372-half*(gtu(1,1)*dd_gtd(1,1,2,2)+gtu(2,2)*dd_
     #gtd(2,2,2,2)+gtu(3,3)*dd_gtd(3,3,2,2)-2*gtd(1,2)*d_Gamt(2,1)-2*gtd
     #(2,2)*d_Gamt(2,2)-2*gtd(2,3)*d_Gamt(2,3))-gtu(1,2)*dd_gtd(1,2,2,2)
     #-gtu(1,3)*dd_gtd(1,3,2,2)-gtu(2,3)*dd_gtd(2,3,2,2)
      Rtd(2,3) = fourth*t448-half*(gtu(1,1)*dd_gtd(1,1,2,3)+gtu(2,2)*dd_
     #gtd(2,2,2,3)+gtu(3,3)*dd_gtd(3,3,2,3)-gtd(1,2)*d_Gamt(3,1)-gtd(1,3
     #)*d_Gamt(2,1)-gtd(2,2)*d_Gamt(3,2)-gtd(2,3)*d_Gamt(2,2)-gtd(2,3)*d
     #_Gamt(3,3)-gtd(3,3)*d_Gamt(2,3))-gtu(1,2)*dd_gtd(1,2,2,3)-gtu(1,3)
     #*dd_gtd(1,3,2,3)-gtu(2,3)*dd_gtd(2,3,2,3)
      Rtd(3,1) = Rtd(1,3)
      Rtd(3,2) = Rtd(2,3)
      Rtd(3,3) = fourth*t527-half*(gtu(1,1)*dd_gtd(1,1,3,3)+gtu(2,2)*dd_
     #gtd(2,2,3,3)+gtu(3,3)*dd_gtd(3,3,3,3)-2*gtd(1,3)*d_Gamt(3,1)-2*gtd
     #(2,3)*d_Gamt(3,2)-2*gtd(3,3)*d_Gamt(3,3))-gtu(1,2)*dd_gtd(1,2,3,3)
     #-gtu(1,3)*dd_gtd(1,3,3,3)-gtu(2,3)*dd_gtd(2,3,3,3)
      t9 = inv_chi*eight_pi_G
      Psi1(1,1) = chi*(Alpha*Rtd(1,1)-dd_Alpha(1,1)+half*(Ct(1,1,1)*d_Al
     #pha(1)+Ct(2,1,1)*d_Alpha(2)+Ct(3,1,1)*d_Alpha(3)))-t9*pTtd_ADM(1,1
     #)+Alpha*Rpd_1(1,1)-d_Alpha(1)*d_chi(1)
      Psi1(1,2) = chi*(Alpha*Rtd(1,2)-dd_Alpha(1,2)+half*(Ct(1,1,2)*d_Al
     #pha(1)+Ct(2,1,2)*d_Alpha(2)+Ct(3,1,2)*d_Alpha(3)))-t9*pTtd_ADM(1,2
     #)+Alpha*Rpd_1(1,2)-half*(d_Alpha(1)*d_chi(2)+d_Alpha(2)*d_chi(1))
      Psi1(1,3) = chi*(Alpha*Rtd(1,3)-dd_Alpha(1,3)+half*(Ct(1,1,3)*d_Al
     #pha(1)+Ct(2,1,3)*d_Alpha(2)+Ct(3,1,3)*d_Alpha(3)))-t9*pTtd_ADM(1,3
     #)+Alpha*Rpd_1(1,3)-half*(d_Alpha(1)*d_chi(3)+d_Alpha(3)*d_chi(1))
      Psi1(2,1) = Psi1(1,2)
      Psi1(2,2) = chi*(Alpha*Rtd(2,2)-dd_Alpha(2,2)+half*(Ct(1,2,2)*d_Al
     #pha(1)+Ct(2,2,2)*d_Alpha(2)+Ct(3,2,2)*d_Alpha(3)))-t9*pTtd_ADM(2,2
     #)+Alpha*Rpd_1(2,2)-d_Alpha(2)*d_chi(2)
      Psi1(2,3) = chi*(Alpha*Rtd(2,3)-dd_Alpha(2,3)+half*(Ct(1,2,3)*d_Al
     #pha(1)+Ct(2,2,3)*d_Alpha(2)+Ct(3,2,3)*d_Alpha(3)))-t9*pTtd_ADM(2,3
     #)+Alpha*Rpd_1(2,3)-half*(d_Alpha(2)*d_chi(3)+d_Alpha(3)*d_chi(2))
      Psi1(3,1) = Psi1(1,3)
      Psi1(3,2) = Psi1(2,3)
      Psi1(3,3) = chi*(Alpha*Rtd(3,3)-dd_Alpha(3,3)+half*(Ct(1,3,3)*d_Al
     #pha(1)+Ct(2,3,3)*d_Alpha(2)+Ct(3,3,3)*d_Alpha(3)))-t9*pTtd_ADM(3,3
     #)+Alpha*Rpd_1(3,3)-d_Alpha(3)*d_chi(3)
      third_trPsi1 = third*(Psi1(1,1)*gtu(1,1)+Psi1(2,2)*gtu(2,2)+Psi1(3
     #,3)*gtu(3,3)+two*(Psi1(1,2)*gtu(1,2)+Psi1(1,3)*gtu(1,3)+Psi1(2,3)*
     #gtu(2,3)))
      Psi1TF(1,1) = Psi1(1,1)-third_trPsi1*gtd(1,1)
      Psi1TF(1,2) = Psi1(1,2)-third_trPsi1*gtd(1,2)
      Psi1TF(1,3) = Psi1(1,3)-third_trPsi1*gtd(1,3)
      Psi1TF(2,1) = Psi1TF(1,2)
      Psi1TF(2,2) = Psi1(2,2)-third_trPsi1*gtd(2,2)
      Psi1TF(2,3) = Psi1(2,3)-third_trPsi1*gtd(2,3)
      Psi1TF(3,1) = Psi1TF(1,3)
      Psi1TF(3,2) = Psi1TF(2,3)
      Psi1TF(3,3) = Psi1(3,3)-third_trPsi1*gtd(3,3)
      gtd_rhs(1,1) = two*(gtd(1,1)*d_Betau(1,1)+gtd(1,2)*d_Betau(1,2)+gt
     #d(1,3)*d_Betau(1,3)-third*gtd(1,1)*div_Beta-Alpha*Atd(1,1))
      gtd_rhs(1,2) = gtd(1,1)*d_Betau(2,1)+gtd(1,2)*d_Betau(1,1)+gtd(1,2
     #)*d_Betau(2,2)+gtd(2,2)*d_Betau(1,2)+gtd(1,3)*d_Betau(2,3)+gtd(2,3
     #)*d_Betau(1,3)-two*(third*gtd(1,2)*div_Beta+Alpha*Atd(1,2))
      gtd_rhs(1,3) = gtd(1,1)*d_Betau(3,1)+gtd(1,3)*d_Betau(1,1)+gtd(1,2
     #)*d_Betau(3,2)+gtd(2,3)*d_Betau(1,2)+gtd(1,3)*d_Betau(3,3)+gtd(3,3
     #)*d_Betau(1,3)-two*(third*gtd(1,3)*div_Beta+Alpha*Atd(1,3))
      gtd_rhs(2,1) = gtd_rhs(1,2)
      gtd_rhs(2,2) = two*(gtd(1,2)*d_Betau(2,1)+gtd(2,2)*d_Betau(2,2)+gt
     #d(2,3)*d_Betau(2,3)-third*gtd(2,2)*div_Beta-Alpha*Atd(2,2))
      gtd_rhs(2,3) = gtd(1,2)*d_Betau(3,1)+gtd(1,3)*d_Betau(2,1)+gtd(2,2
     #)*d_Betau(3,2)+gtd(2,3)*d_Betau(2,2)+gtd(2,3)*d_Betau(3,3)+gtd(3,3
     #)*d_Betau(2,3)-two*(third*gtd(2,3)*div_Beta+Alpha*Atd(2,3))
      gtd_rhs(3,1) = gtd_rhs(1,3)
      gtd_rhs(3,2) = gtd_rhs(2,3)
      gtd_rhs(3,3) = two*(gtd(1,3)*d_Betau(3,1)+gtd(2,3)*d_Betau(3,2)+gt
     #d(3,3)*d_Betau(3,3)-third*gtd(3,3)*div_Beta-Alpha*Atd(3,3))
      t14 = Alpha*trK-twothirds*div_Beta
      t24 = two*Alpha
      Atd_rhs(1,1) = two*(Atd(1,1)*(d_Betau(1,1)-Alpha*Atud(1,1))+Atd(1,
     #2)*(d_Betau(1,2)-Alpha*Atud(2,1))+Atd(1,3)*(d_Betau(1,3)-Alpha*Atu
     #d(3,1)))+Atd(1,1)*t14+Psi1TF(1,1)
      Atd_rhs(1,2) = Atd(1,1)*d_Betau(2,1)+Atd(1,2)*d_Betau(1,1)+Atd(1,2
     #)*d_Betau(2,2)+Atd(2,2)*d_Betau(1,2)+Atd(1,3)*d_Betau(2,3)+Atd(2,3
     #)*d_Betau(1,3)+Atd(1,2)*t14+Psi1TF(1,2)-t24*(Atd(1,1)*Atud(1,2)+At
     #d(1,2)*Atud(2,2)+Atd(1,3)*Atud(3,2))
      Atd_rhs(1,3) = Atd(1,1)*d_Betau(3,1)+Atd(1,3)*d_Betau(1,1)+Atd(1,2
     #)*d_Betau(3,2)+Atd(2,3)*d_Betau(1,2)+Atd(1,3)*d_Betau(3,3)+Atd(3,3
     #)*d_Betau(1,3)+Atd(1,3)*t14+Psi1TF(1,3)-t24*(Atd(1,1)*Atud(1,3)+At
     #d(1,2)*Atud(2,3)+Atd(1,3)*Atud(3,3))
      Atd_rhs(2,1) = Atd_rhs(1,2)
      Atd_rhs(2,2) = two*(Atd(1,2)*(d_Betau(2,1)-Alpha*Atud(1,2))+Atd(2,
     #2)*(d_Betau(2,2)-Alpha*Atud(2,2))+Atd(2,3)*(d_Betau(2,3)-Alpha*Atu
     #d(3,2)))+Atd(2,2)*t14+Psi1TF(2,2)
      Atd_rhs(2,3) = Atd(1,2)*d_Betau(3,1)+Atd(1,3)*d_Betau(2,1)+Atd(2,2
     #)*d_Betau(3,2)+Atd(2,3)*d_Betau(2,2)+Atd(2,3)*d_Betau(3,3)+Atd(3,3
     #)*d_Betau(2,3)+Atd(2,3)*t14+Psi1TF(2,3)-t24*(Atd(1,2)*Atud(1,3)+At
     #d(2,2)*Atud(2,3)+Atd(2,3)*Atud(3,3))
      Atd_rhs(3,1) = Atd_rhs(1,3)
      Atd_rhs(3,2) = Atd_rhs(2,3)
      Atd_rhs(3,3) = two*(Atd(1,3)*(d_Betau(3,1)-Alpha*Atud(1,3))+Atd(2,
     #3)*(d_Betau(3,2)-Alpha*Atud(2,3))+Atd(3,3)*(d_Betau(3,3)-Alpha*Atu
     #d(3,3)))+Atd(3,3)*t14+Psi1TF(3,3)
      t28 = eight_pi_G*inv_chi
      t30 = twothirds*d_trK(1)+t28*Jtd_ADM(1)
      t34 = twothirds*d_trK(2)+t28*Jtd_ADM(2)
      t38 = twothirds*d_trK(3)+t28*Jtd_ADM(3)
      t42 = threehalves*Alpha
      t45 = -t42*inv_chi*d_chi(1)-d_Alpha(1)
      t49 = -t42*inv_chi*d_chi(2)-d_Alpha(2)
      t53 = -t42*inv_chi*d_chi(3)-d_Alpha(3)
      Gamt_rhs(1) = twothirds*CalGamt(1)*div_Beta-CalGamt(1)*d_Betau(1,1
     #)-CalGamt(2)*d_Betau(2,1)-CalGamt(3)*d_Betau(3,1)+gtu(1,1)*dd_Beta
     #u(1,1,1)+gtu(2,2)*dd_Betau(2,2,1)+gtu(3,3)*dd_Betau(3,3,1)+two*(gt
     #u(1,2)*dd_Betau(1,2,1)+gtu(1,3)*dd_Betau(1,3,1)+gtu(2,3)*dd_Betau(
     #2,3,1))+third*(gtu(1,1)*d_div_Beta(1)+gtu(1,2)*d_div_Beta(2)+gtu(1
     #,3)*d_div_Beta(3))+two*(Alpha*(half*(Ct(1,1,1)*Atu(1,1)+Ct(1,2,2)*
     #Atu(2,2)+Ct(1,3,3)*Atu(3,3))+Ct(1,1,2)*Atu(1,2)+Ct(1,1,3)*Atu(1,3)
     #+Ct(1,2,3)*Atu(2,3)-gtu(1,1)*t30-gtu(1,2)*t34-gtu(1,3)*t38)+Atu(1,
     #1)*t45+Atu(1,2)*t49+Atu(1,3)*t53)
      Gamt_rhs(2) = twothirds*CalGamt(2)*div_Beta-CalGamt(1)*d_Betau(1,2
     #)-CalGamt(2)*d_Betau(2,2)-CalGamt(3)*d_Betau(3,2)+gtu(1,1)*dd_Beta
     #u(1,1,2)+gtu(2,2)*dd_Betau(2,2,2)+gtu(3,3)*dd_Betau(3,3,2)+two*(gt
     #u(1,2)*dd_Betau(1,2,2)+gtu(1,3)*dd_Betau(1,3,2)+gtu(2,3)*dd_Betau(
     #2,3,2))+third*(gtu(1,2)*d_div_Beta(1)+gtu(2,2)*d_div_Beta(2)+gtu(2
     #,3)*d_div_Beta(3))+two*(Alpha*(half*(Ct(2,1,1)*Atu(1,1)+Ct(2,2,2)*
     #Atu(2,2)+Ct(2,3,3)*Atu(3,3))+Ct(2,1,2)*Atu(1,2)+Ct(2,1,3)*Atu(1,3)
     #+Ct(2,2,3)*Atu(2,3)-gtu(1,2)*t30-gtu(2,2)*t34-gtu(2,3)*t38)+Atu(1,
     #2)*t45+Atu(2,2)*t49+Atu(2,3)*t53)
      Gamt_rhs(3) = twothirds*CalGamt(3)*div_Beta-CalGamt(1)*d_Betau(1,3
     #)-CalGamt(2)*d_Betau(2,3)-CalGamt(3)*d_Betau(3,3)+gtu(1,1)*dd_Beta
     #u(1,1,3)+gtu(2,2)*dd_Betau(2,2,3)+gtu(3,3)*dd_Betau(3,3,3)+two*(gt
     #u(1,2)*dd_Betau(1,2,3)+gtu(1,3)*dd_Betau(1,3,3)+gtu(2,3)*dd_Betau(
     #2,3,3))+third*(gtu(1,3)*d_div_Beta(1)+gtu(2,3)*d_div_Beta(2)+gtu(3
     #,3)*d_div_Beta(3))+two*(Alpha*(half*(Ct(3,1,1)*Atu(1,1)+Ct(3,2,2)*
     #Atu(2,2)+Ct(3,3,3)*Atu(3,3))+Ct(3,1,2)*Atu(1,2)+Ct(3,1,3)*Atu(1,3)
     #+Ct(3,2,3)*Atu(2,3)-gtu(1,3)*t30-gtu(2,3)*t34-gtu(3,3)*t38)+Atu(1,
     #3)*t45+Atu(2,3)*t49+Atu(3,3)*t53)
      t3 = threefourths*(lambda_f0+lambda_f1*Alpha)
      Betau_rhs(1) = t3*Bu(1)
      Betau_rhs(2) = t3*Bu(2)
      Betau_rhs(3) = t3*Bu(3)
      Bu_rhs(1) = Gamt_rhs(1)-Bu(1)*feta
      Bu_rhs(2) = Gamt_rhs(2)-Bu(2)*feta
      Bu_rhs(3) = Gamt_rhs(3)-Bu(3)*feta
      Alpha_rhs = -two*Alpha*(trK-trK0)
      chi_rhs = -twothirds*chi*(div_Beta-Alpha*trK)
      t28 = trK**2
      t34 = gtu(1,1)
      t37 = gtu(1,2)
      t41 = gtu(1,3)
      t45 = gtu(2,2)
      t48 = gtu(2,3)
      t52 = gtu(3,3)
      t56 = d_Alpha(1)
      t59 = d_Alpha(2)
      t62 = d_Alpha(3)
      t67 = d_chi(1)
      t70 = d_chi(2)
      t73 = d_chi(3)

      trK_rhs = Alpha*(Atd(1,1)*Atu(1,1)+Atd(2,2)*Atu(2,2)+Atd(3,3)*Atu(
     #3,3)+two*(Atd(1,2)*Atu(1,2)+Atd(1,3)*Atu(1,3)+Atd(2,3)*Atu(2,3))+t
     #hird*t28+four_pi_G*(rho_ADM+tr_pT))-chi*(t34*dd_Alpha(1,1)+2*t37*d
     #d_Alpha(1,2)+2*t41*dd_Alpha(1,3)+t45*dd_Alpha(2,2)+2*t48*dd_Alpha(
     #2,3)+t52*dd_Alpha(3,3)-CalGamt(1)*t56-CalGamt(2)*t59-CalGamt(3)*t6
     #2)+half*(t34*t56*t67+t37*t56*t70+t41*t56*t73+t37*t59*t67+t45*t59*t
     #70+t48*t59*t73+t41*t62*t67+t48*t62*t70+t52*t62*t73)

      trK2_rhs = Alpha*(Atd(1,1)*Atu(1,1)+Atd(2,2)*Atu(2,2)+Atd(3,3)*Atu(
     #3,3)+two*(Atd(1,2)*Atu(1,2)+Atd(1,3)*Atu(1,3)+Atd(2,3)*Atu(2,3))+t
     #hird*t28+four_pi_G*(rho_ADM+tr_pT))

      trK3_rhs = -chi*(t34*dd_Alpha(1,1)+2*t37*d
     #d_Alpha(1,2)+2*t41*dd_Alpha(1,3)+t45*dd_Alpha(2,2)+2*t48*dd_Alpha(
     #2,3)+t52*dd_Alpha(3,3)-CalGamt(1)*t56-CalGamt(2)*t59-CalGamt(3)*t6
     #2)+half*(t34*t56*t67+t37*t56*t70+t41*t56*t73+t37*t59*t67+t45*t59*t
     #70+t48*t59*t73+t41*t62*t67+t48*t62*t70+t52*t62*t73)

#if 0
      t5 = d_chi(1)
      t8 = d_chi(2)
      t11 = d_chi(3)
      divE = d_emEu(1,1)+d_emEu(2,2)+d_emEu(3,3)-threehalves*(emEu(1)*t5
     #+emEu(2)*t8+emEu(3)*t11)*inv_chi
      t19 = emBu(1)
      t21 = emBu(2)
      t23 = emBu(3)
      divB = d_emBu(1,1)+d_emBu(2,2)+d_emBu(3,3)-threehalves*(t19*t5+t21
     #*t8+t23*t11)*inv_chi
      t28 = t19**2
      t31 = t21**2
      t34 = t23**2
      Bsq = t28*gtd(1,1)+t31*gtd(2,2)+t34*gtd(3,3)+two*(t19*t21*gtd(1,2)
     #+t19*t23*gtd(1,3)+t21*t23*gtd(2,3))
      inv_Bsq = one/Bsq
      t1 = divE*sqrt_chi
      t4 = -emEu(3)*emBu(2)+emEu(2)*emBu(3)
      t8 = emEu(3)*emBu(1)-emEu(1)*emBu(3)
      t12 = -emEu(2)*emBu(1)+emEu(1)*emBu(2)
      Ju(1) = t1*(gtu(1,1)*t4+gtu(1,2)*t8+gtu(1,3)*t12)*inv_Bsq
      Ju(2) = t1*(gtu(1,2)*t4+gtu(2,2)*t8+gtu(2,3)*t12)*inv_Bsq
      Ju(3) = t1*(gtu(1,3)*t4+gtu(2,3)*t8+gtu(3,3)*t12)*inv_Bsq
      !set by hand to zero
      Ju(1) = 0.0
      Ju(2) = 0.0
      Ju(3) = 0.0
      t7 = emBu(1)*gtd(1,2)+emBu(2)*gtd(2,2)+emBu(3)*gtd(2,3)
      t10 = d_Alpha(3)-Alpha*d_chi(3)*inv_chi
      t23 = emBu(1)*gtd(1,3)+emBu(2)*gtd(2,3)+emBu(3)*gtd(3,3)
      t26 = d_Alpha(2)-Alpha*d_chi(2)*inv_chi
      t53 = emBu(1)*gtd(1,1)+emBu(2)*gtd(1,2)+emBu(3)*gtd(1,3)
      t65 = d_Alpha(1)-Alpha*d_chi(1)*inv_chi
      emEu_rhs(1) = -emEu(1)*d_Betau(1,1)-emEu(2)*d_Betau(2,1)-emEu(3)*d
     #_Betau(3,1)+sqrt_chi*(-t7*t10-Alpha*(d_emBu(3,1)*gtd(1,2)+emBu(1)*
     #d_gtd(3,1,2)+d_emBu(3,2)*gtd(2,2)+emBu(2)*d_gtd(3,2,2)+d_emBu(3,3)
     #*gtd(2,3)+emBu(3)*d_gtd(3,2,3))+t23*t26+Alpha*(d_emBu(2,1)*gtd(1,3
     #)+emBu(1)*d_gtd(2,1,3)+d_emBu(2,2)*gtd(2,3)+emBu(2)*d_gtd(2,2,3)+d
     #_emBu(2,3)*gtd(3,3)+emBu(3)*d_gtd(2,3,3)))+Alpha*(trK*emEu(1)-chi*
     #(gtu(1,1)*d_emPsi(1)+gtu(1,2)*d_emPsi(2)+gtu(1,3)*d_emPsi(3))-Ju(1
     #))
      emEu_rhs(2) = -emEu(1)*d_Betau(1,2)-emEu(2)*d_Betau(2,2)-emEu(3)*d
     #_Betau(3,2)+sqrt_chi*(t53*t10+Alpha*(d_emBu(3,1)*gtd(1,1)+emBu(1)*
     #d_gtd(3,1,1)+d_emBu(3,2)*gtd(1,2)+emBu(2)*d_gtd(3,1,2)+d_emBu(3,3)
     #*gtd(1,3)+emBu(3)*d_gtd(3,1,3))-t23*t65-Alpha*(d_emBu(1,1)*gtd(1,3
     #)+emBu(1)*d_gtd(1,1,3)+d_emBu(1,2)*gtd(2,3)+emBu(2)*d_gtd(1,2,3)+d
     #_emBu(1,3)*gtd(3,3)+emBu(3)*d_gtd(1,3,3)))+Alpha*(trK*emEu(2)-chi*
     #(gtu(1,2)*d_emPsi(1)+gtu(2,2)*d_emPsi(2)+gtu(2,3)*d_emPsi(3))-Ju(2
     #))
      emEu_rhs(3) = -emEu(1)*d_Betau(1,3)-emEu(2)*d_Betau(2,3)-emEu(3)*d
     #_Betau(3,3)+sqrt_chi*(-t53*t26-Alpha*(d_emBu(2,1)*gtd(1,1)+emBu(1)
     #*d_gtd(2,1,1)+d_emBu(2,2)*gtd(1,2)+emBu(2)*d_gtd(2,1,2)+d_emBu(2,3
     #)*gtd(1,3)+emBu(3)*d_gtd(2,1,3))+t7*t65+Alpha*(d_emBu(1,1)*gtd(1,2
     #)+emBu(1)*d_gtd(1,1,2)+d_emBu(1,2)*gtd(2,2)+emBu(2)*d_gtd(1,2,2)+d
     #_emBu(1,3)*gtd(2,3)+emBu(3)*d_gtd(1,2,3)))+Alpha*(trK*emEu(3)-chi*
     #(gtu(1,3)*d_emPsi(1)+gtu(2,3)*d_emPsi(2)+gtu(3,3)*d_emPsi(3))-Ju(3
     #))
      t7 = emEu(1)*gtd(1,2)+emEu(2)*gtd(2,2)+emEu(3)*gtd(2,3)
      t10 = d_Alpha(3)-Alpha*d_chi(3)*inv_chi
      t23 = emEu(1)*gtd(1,3)+emEu(2)*gtd(2,3)+emEu(3)*gtd(3,3)
      t26 = d_Alpha(2)-Alpha*d_chi(2)*inv_chi
      t53 = emEu(1)*gtd(1,1)+emEu(2)*gtd(1,2)+emEu(3)*gtd(1,3)
      t65 = d_Alpha(1)-Alpha*d_chi(1)*inv_chi
      emBu_rhs(1) = -emBu(1)*d_Betau(1,1)-emBu(2)*d_Betau(2,1)-emBu(3)*d
     #_Betau(3,1)-sqrt_chi*(-t7*t10-Alpha*(d_emEu(3,1)*gtd(1,2)+emEu(1)*
     #d_gtd(3,1,2)+d_emEu(3,2)*gtd(2,2)+emEu(2)*d_gtd(3,2,2)+d_emEu(3,3)
     #*gtd(2,3)+emEu(3)*d_gtd(3,2,3))+t23*t26+Alpha*(d_emEu(2,1)*gtd(1,3
     #)+emEu(1)*d_gtd(2,1,3)+d_emEu(2,2)*gtd(2,3)+emEu(2)*d_gtd(2,2,3)+d
     #_emEu(2,3)*gtd(3,3)+emEu(3)*d_gtd(2,3,3)))+Alpha*(trK*emBu(1)+chi*
     #(gtu(1,1)*d_emPhi(1)+gtu(1,2)*d_emPhi(2)+gtu(1,3)*d_emPhi(3)))
      emBu_rhs(2) = -emBu(1)*d_Betau(1,2)-emBu(2)*d_Betau(2,2)-emBu(3)*d
     #_Betau(3,2)-sqrt_chi*(t53*t10+Alpha*(d_emEu(3,1)*gtd(1,1)+emEu(1)*
     #d_gtd(3,1,1)+d_emEu(3,2)*gtd(1,2)+emEu(2)*d_gtd(3,1,2)+d_emEu(3,3)
     #*gtd(1,3)+emEu(3)*d_gtd(3,1,3))-t23*t65-Alpha*(d_emEu(1,1)*gtd(1,3
     #)+emEu(1)*d_gtd(1,1,3)+d_emEu(1,2)*gtd(2,3)+emEu(2)*d_gtd(1,2,3)+d
     #_emEu(1,3)*gtd(3,3)+emEu(3)*d_gtd(1,3,3)))+Alpha*(trK*emBu(2)+chi*
     #(gtu(1,2)*d_emPhi(1)+gtu(2,2)*d_emPhi(2)+gtu(2,3)*d_emPhi(3)))
      emBu_rhs(3) = -emBu(1)*d_Betau(1,3)-emBu(2)*d_Betau(2,3)-emBu(3)*d
     #_Betau(3,3)-sqrt_chi*(-t53*t26-Alpha*(d_emEu(2,1)*gtd(1,1)+emEu(1)
     #*d_gtd(2,1,1)+d_emEu(2,2)*gtd(1,2)+emEu(2)*d_gtd(2,1,2)+d_emEu(2,3
     #)*gtd(1,3)+emEu(3)*d_gtd(2,1,3))+t7*t65+Alpha*(d_emEu(1,1)*gtd(1,2
     #)+emEu(1)*d_gtd(1,1,2)+d_emEu(1,2)*gtd(2,2)+emEu(2)*d_gtd(1,2,2)+d
     #_emEu(1,3)*gtd(2,3)+emEu(3)*d_gtd(1,2,3)))+Alpha*(trK*emBu(3)+chi*
     #(gtu(1,3)*d_emPhi(1)+gtu(2,3)*d_emPhi(2)+gtu(3,3)*d_emPhi(3)))
      emPsi_rhs = -kappa_1*Alpha*emPsi
      emPhi_rhs = Alpha*(divB-kappa_2*emPhi)
#endif
