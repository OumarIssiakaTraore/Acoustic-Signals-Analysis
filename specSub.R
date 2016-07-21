##' Noise cancelation by spectral substraction
##' 
##' 
##' @encoding latin1
##' @title Noise cancelation by spectral substraction
##'  
##' @param signal detected raw signal
##' @param n number of subdivision to be considered (to be define considering the sampling frequency between 100-1000 for 20kHz)
##' @param f sampling frequency
##' @param alpha over-substraction factor (1 means no over-substraction)
##' @return Denoised signal
##' 
##' 
##' @note For instant, this fuction is specific to our research topic. It have to be adaptated.
##' @author TRAORE Oumar Issiaka
##' 
##' 
##' 
##' 
##' 
##' 
##' 



##-------------------------------------------------------------------------------------------------
## Creation date: Octobre 2015
## Last modification date: 19/01/2016
##-------------------------------------------------------------------------------------------------


  ####----- Function start here

    specSub<-function(signal,time=NULL,n,f,alpha,plot=FALSE){
  

  
  ##
  ##--  Creating the data frame for the processing
  ##
      
      leng<-length(signal)
      if(is.null(time)){
        time<-as.data.frame(seq(0,leng,length.out=leng)/f)
      }
      data<-as.data.frame(cbind(time,signal))
      colnames(data)[1]<-'time'

  
  ##  
  ##--- Dividing the signal in subsignals
  ##
  
    ind<- cut(round(1:length(signal)), quantile(round(1:length(signal)), probs = seq(0, 1, by = (1/n))))
    subSignals<-split(data,ind)
    
   

  ##
  ##-- Correction vectors lengths
  ##

    lCorrection<-min(do.call(rbind,lapply(subSignals,nrow)))
                    
    subSignals<-lapply(subSignals,FUN= function(x,lCorrection){
      return(x[1:lCorrection,])
             },lCorrection)  
  
    
  
  ##
  ##-- Windowing sub-signals and computing fast fourier transforms 
  ##

    fftsubSignals<-lapply(subSignals,FUN=function(x){
      w <- ftwindow(nrow(x), wn = "hanning")
      TF<-fft(x$signal*w)/nrow(x)
      return(TF)
      
    })


  ##
  ##---  Spectral substraction and inverse fourier transform of sub-signals
  ##

    noise<-fftsubSignals[[1]]
    #x<-fftsubSignals[[4]]
    
    FinalSubSignals<-lapply(fftsubSignals,FUN=function(x,alpha,noise){
      
      Y<-vector()
      phase<-Arg(x)
      modS<-Mod(x)
      modN<-alpha*Mod(noise)
      
    Y<- ifelse((modS-modN)<0,0,(modS-modN)*exp(phase*complex(imaginary=1)) )
    Z<-fft(Y,inv=TRUE)
    Z<-as.data.frame(Re(Z)) 
    colnames(Z)<-c('Denoised_Signal')
      return(Z)
    
      
    },alpha,noise=fftsubSignals[[1]])

    
  
  ##
  ##--  Signal reconstruction
  ##
  
    FinalSubSignals<-do.call(rbind,FinalSubSignals)
    subSignals<-do.call(rbind,subSignals)
    resultat<-cbind(subSignals,FinalSubSignals)

  





  ##
  ##--  Ploting the processing result
  ##
    if(plot==TRUE){
    plot(resultat$time,resultat$signal,type='l',col='#798081',xlab='time(s)',ylab='Amplitude (V)',
         main='Waveform comparison')
    lines(resultat$time,resultat$Denoised_Signal,col="#CDAA7D")
    legend('topright',legend=c('Noised Signal','Denoised Signal'),col=c('#798081',"#CDAA7D"),pch=15)
    legend('bottomright',c(paste('alpha = ',alpha),
                           paste("n = ",n),
                           paste('sampling frequency = ',f)))
    grid(20,20)
    }
  ##
  ##-- Return the processing result
  ##
    
    return(resultat)



}




