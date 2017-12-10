      Gamt_rhs(1) = Gamt_rhs(1)+Betau(1)*adv_d_Gamt(1,1)+Betau(2)*ad
     &v_d_Gamt(2,1)+Betau(3)*adv_d_Gamt(3,1)
      Gamt_rhs(2) = Gamt_rhs(2)+Betau(1)*adv_d_Gamt(1,2)+Betau(2)*ad
     &v_d_Gamt(2,2)+Betau(3)*adv_d_Gamt(3,2)
      Gamt_rhs(3) = Gamt_rhs(3)+Betau(1)*adv_d_Gamt(1,3)+Betau(2)*ad
     &v_d_Gamt(2,3)+Betau(3)*adv_d_Gamt(3,3)
      Atd_rhs(1,1) = Atd_rhs(1,1)+Betau(1)*adv_d_Atd(1,1,1)+Betau(2)
     &*adv_d_Atd(2,1,1)+Betau(3)*adv_d_Atd(3,1,1)
      Atd_rhs(1,2) = Atd_rhs(1,2)+Betau(1)*adv_d_Atd(1,1,2)+Betau(2)
     &*adv_d_Atd(2,1,2)+Betau(3)*adv_d_Atd(3,1,2)
      Atd_rhs(1,3) = Atd_rhs(1,3)+Betau(1)*adv_d_Atd(1,1,3)+Betau(2)
     &*adv_d_Atd(2,1,3)+Betau(3)*adv_d_Atd(3,1,3)
      Atd_rhs(2,1) = Atd_rhs(1,2)
      Atd_rhs(2,2) = Atd_rhs(2,2)+Betau(1)*adv_d_Atd(1,2,2)+Betau(2)
     &*adv_d_Atd(2,2,2)+Betau(3)*adv_d_Atd(3,2,2)
      Atd_rhs(2,3) = Atd_rhs(2,3)+Betau(1)*adv_d_Atd(1,2,3)+Betau(2)
     &*adv_d_Atd(2,2,3)+Betau(3)*adv_d_Atd(3,2,3)
      Atd_rhs(3,1) = Atd_rhs(1,3)
      Atd_rhs(3,2) = Atd_rhs(2,3)
      Atd_rhs(3,3) = Atd_rhs(3,3)+Betau(1)*adv_d_Atd(1,3,3)+Betau(2)
     &*adv_d_Atd(2,3,3)+Betau(3)*adv_d_Atd(3,3,3)
      gtd_rhs(1,1) = gtd_rhs(1,1)+Betau(1)*adv_d_gtd(1,1,1)+Betau(2)
     &*adv_d_gtd(2,1,1)+Betau(3)*adv_d_gtd(3,1,1)
      gtd_rhs(1,2) = gtd_rhs(1,2)+Betau(1)*adv_d_gtd(1,1,2)+Betau(2)
     &*adv_d_gtd(2,1,2)+Betau(3)*adv_d_gtd(3,1,2)
      gtd_rhs(1,3) = gtd_rhs(1,3)+Betau(1)*adv_d_gtd(1,1,3)+Betau(2)
     &*adv_d_gtd(2,1,3)+Betau(3)*adv_d_gtd(3,1,3)
      gtd_rhs(2,1) = gtd_rhs(1,2)
      gtd_rhs(2,2) = gtd_rhs(2,2)+Betau(1)*adv_d_gtd(1,2,2)+Betau(2)
     &*adv_d_gtd(2,2,2)+Betau(3)*adv_d_gtd(3,2,2)
      gtd_rhs(2,3) = gtd_rhs(2,3)+Betau(1)*adv_d_gtd(1,2,3)+Betau(2)
     &*adv_d_gtd(2,2,3)+Betau(3)*adv_d_gtd(3,2,3)
      gtd_rhs(3,1) = gtd_rhs(1,3)
      gtd_rhs(3,2) = gtd_rhs(2,3)
      gtd_rhs(3,3) = gtd_rhs(3,3)+Betau(1)*adv_d_gtd(1,3,3)+Betau(2)
     &*adv_d_gtd(2,3,3)+Betau(3)*adv_d_gtd(3,3,3)
      Betau_rhs(1) = Betau_rhs(1)+lambda_2*(Betau(1)*adv_d_Betau(1,1
     &)+Betau(2)*adv_d_Betau(2,1)+Betau(3)*adv_d_Betau(3,1))
      Betau_rhs(2) = Betau_rhs(2)+lambda_2*(Betau(1)*adv_d_Betau(1,2
     &)+Betau(2)*adv_d_Betau(2,2)+Betau(3)*adv_d_Betau(3,2))
      Betau_rhs(3) = Betau_rhs(3)+lambda_2*(Betau(1)*adv_d_Betau(1,3
     &)+Betau(2)*adv_d_Betau(2,3)+Betau(3)*adv_d_Betau(3,3))
      t2 = 1-lambda_4
      Bu_rhs(1) = Bu_rhs(1)+Betau(1)*(lambda_3*adv_d_Bu(1,1)+t2*adv_
     &d_Gamt(1,1))+Betau(2)*(lambda_3*adv_d_Bu(2,1)+t2*adv_d_Gamt(2,1))+
     &Betau(3)*(lambda_3*adv_d_Bu(3,1)+t2*adv_d_Gamt(3,1))
      Bu_rhs(2) = Bu_rhs(2)+Betau(1)*(lambda_3*adv_d_Bu(1,2)+t2*adv_
     &d_Gamt(1,2))+Betau(2)*(lambda_3*adv_d_Bu(2,2)+t2*adv_d_Gamt(2,2))+
     &Betau(3)*(lambda_3*adv_d_Bu(3,2)+t2*adv_d_Gamt(3,2))
      Bu_rhs(3) = Bu_rhs(3)+Betau(1)*(lambda_3*adv_d_Bu(1,3)+t2*adv_
     &d_Gamt(1,3))+Betau(2)*(lambda_3*adv_d_Bu(2,3)+t2*adv_d_Gamt(2,3))+
     &Betau(3)*(lambda_3*adv_d_Bu(3,3)+t2*adv_d_Gamt(3,3))
      t1 = Betau(1)
      t4 = Betau(2)
      t7 = Betau(3)
      chi_rhs = chi_rhs+t1*adv_d_chi(1)+t4*adv_d_chi(2)+t7*adv_d_chi
     &(3)
      trK_rhs = trK_rhs+t1*adv_d_trK(1)+t4*adv_d_trK(2)+t7*adv_d_trK
     &(3)
      Alpha_rhs = Alpha_rhs+lambda_1*(t1*adv_d_Alpha(1)+t4*adv_d_Alp
     &ha(2)+t7*adv_d_Alpha(3))

#if 0
      emEu_rhs(1) = emEu_rhs(1)+Betau(1)*adv_d_emEu(1,1)+Betau(2)*ad
     &v_d_emEu(2,1)+Betau(3)*adv_d_emEu(3,1)
      emEu_rhs(2) = emEu_rhs(2)+Betau(1)*adv_d_emEu(1,2)+Betau(2)*ad
     &v_d_emEu(2,2)+Betau(3)*adv_d_emEu(3,2)
      emEu_rhs(3) = emEu_rhs(3)+Betau(1)*adv_d_emEu(1,3)+Betau(2)*ad
     &v_d_emEu(2,3)+Betau(3)*adv_d_emEu(3,3)
      emBu_rhs(1) = emBu_rhs(1)+Betau(1)*adv_d_emBu(1,1)+Betau(2)*ad
     &v_d_emBu(2,1)+Betau(3)*adv_d_emBu(3,1)
      emBu_rhs(2) = emBu_rhs(2)+Betau(1)*adv_d_emBu(1,2)+Betau(2)*ad
     &v_d_emBu(2,2)+Betau(3)*adv_d_emBu(3,2)
      emBu_rhs(3) = emBu_rhs(3)+Betau(1)*adv_d_emBu(1,3)+Betau(2)*ad
     &v_d_emBu(2,3)+Betau(3)*adv_d_emBu(3,3)
      t1 = Betau(1)
      t4 = Betau(2)
      t7 = Betau(3)
      emPsi_rhs = emPsi_rhs+t1*adv_d_emPsi(1)+t4*adv_d_emPsi(2)+t7*a
     &dv_d_emPsi(3)
      emPhi_rhs = emPhi_rhs+t1*adv_d_emPhi(1)+t4*adv_d_emPhi(2)+t7*a
     &dv_d_emPhi(3)
#endif
