

Turc<-function(A,B,C,D,E,F,
                    NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  v<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  Actual_Eva<-(prec)*(((1+(Pot_Eva/(prec))^(-v)))^(-1/v))
  
  
  return(Actual_Eva)
}




Fu<-function(A,B,C,D,E,F,
             NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  
  omega<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  Actual_Eva<-(prec)*(1+(Pot_Eva/(prec))-(1+(Pot_Eva/(prec))^(omega))^(1/omega))
  
  return(Actual_Eva)
}

Schreiber<-function(prec,Pot_Eva)
{
  Actual_Eva<-(prec)*(1-exp(-1*(Pot_Eva/(prec))))
  return(Actual_Eva)
}

Oldekop<-function(prec,Pot_Eva)
{
  Actual_Eva<-(prec)*(Pot_Eva/prec)*tanh(1/(Pot_Eva/prec))
  return(Actual_Eva)
}

Budyko<-function(prec,Pot_Eva)
{
  phi<-(Pot_Eva/prec)
  
  a<-0.5
  Actual_Eva<-(prec)*((phi*tanh(1/phi)*(1-exp(-1*phi)))^a)
  
  return(Actual_Eva)
}

Miley_alpha<-function(A,B,C,D,E,F,
                      NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  alpha<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-(prec)*(exp(alpha*(1-1/phi))-1)/(exp(alpha*(1-1/phi))-1/phi)
  
  return(Actual_Eva)
}


Miley_beta<-function(A,B,C,D,E,F,
                      NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  beta<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-(prec)*(exp(beta*(phi-1))-1)/(exp(beta*(phi-1))-1/phi)
  

  return(Actual_Eva)
}


Milly<-function(A,B,C,D,E,F,
                NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  W<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-(prec)*((1+W*phi)/(1+W*phi+1/phi))
  
  
  return(Actual_Eva)
}

Zhang<-function(A,B,C,D,E,F,
                NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  W<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-(prec)*((1+phi)-((1+(phi^W))^(1/W)))
  
  return(Actual_Eva)
}




Linear<-function(A,B,C,D,E,F,
                 NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  b<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-b*phi
  
  return(Actual_Eva)
}


Zhou<-function(A,B,C,D,E,F,
               n,
               NDVI,Long,Lat,Alt,curve_number,prec,Pot_Eva)
{
  k<-(A*NDVI)+(B*Long)+(C*Lat)+(D*Alt)+(E*curve_number)+F
  
  
  
  phi<-(Pot_Eva/prec)
  
  Actual_Eva<-phi/((1+k*(phi^n))^(1/n))
  
  return(Actual_Eva)
}
