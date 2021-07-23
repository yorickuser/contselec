# file contselec/R/contselec-package.R
# copyright (C) 2018 Hiroshi C. Ito
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 2 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/



#' Model selection procesure by contracting explanatry variables of high absolute correlation, with examining statistical strength of each explanatry variable.
#' @aliases contselec contselec-package
#' @importFrom foreach %do%
#' @keywords internal
"_PACKAGE"

flag_envs=FALSE;
flag_X11=FALSE;

#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/ FUNCTION FOR CONTRACTION AND MODEL SELECTION  _/_/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
##' This function contracts explanatory variables and conduct best subset model selection, permutation test and stepwise model selection.
##'
##' Spatial correlation of residuals of each fitted model can also be examimed.
##' @title contracts explanatory variables, conducts model selection and examins statistical strength of explanatory variables
##' @param data00 data.frame : an explained variable and explanatory variables. The explained variable is specified by "target".
##' @param pos_x array of numeric : x-coordinates of spatial position for each sampled data in "data00".
##' @param pos_y array of numeric : y-coordinates of spatial position for each sampled data in "data00".
##' @param target integer or character : spefifies explained variables by its column id (integer) or column name (character).
##' @param edge_cor real number : threshold for grouping; variables of correlation no less than this threshold are grouped together.
##' @param edge_explain real number : threshold for the variance ratio to be explained, based on which the number of principal components "np" to represent the contracted group is determined.
##' @param edge_param_number real number: this parameter constrains the number of explanatory variables in the subset model selection and stepwise model selection, so that [sample size]/[number of free parameter] >= "edge_param_number".
##' @param repeatn integer : the number of model evaluations proccessed as a single job for each CPU core (This parameter is meaningful only when "use_pforeach" is TRUE).
##' @param family character : error distribution. Possible choices are "glmm_poisson", "poisson","gaussian".
##' @param check_spa_cor TRUE/FALSE : if TRUE, the best model with residuals significantly correlated with c(pos_x, pos_y) are removed repeatedly until the correlation becomes non-significant.
##' @param use_pforeach  TRUE/FALSE : if TRUE, pforeach is used instead of foreach in best subset selection and permutation tests (when "perm" is TRUE).
##' @param perm  TRUE/FALSE : if TRUE, permutation tests for explained variables are conducted.
##' @param nrep  integer : number of resampling in the permutation test.
##' @return list(target=target, result=result, data0=data0, data1=data1, group=group, family=family, bestmodels=res_sb$bestmodels, bestmodel_stepwise=res_ss$bestmodel, bestmodel=res_sb$bestmodel, pos_x=pos_x, pos_y=pos_y)
##' @author Hiroshi C. Ito
##' @examples
##' 
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' 
##' res=cont_selec(data,target="Horsepower",edge_cor=0.9,edge_explain=0.6,
##'                edge_param_number=3,family="gaussian",use_pforeach=FALSE);
##'
##' plot_each_effect(res)
##' plot_combined_effect(res,sign_effect=-1)
##' plot_combined_effect(res,sign_effect=1)
##' 
##' print(res);
##' summary(res);
#' @export
cont_selec <- function(data00,pos_x=NULL,pos_y=NULL,target=1,edge_cor=0.52,edge_explain=0.65,edge_param_number=3,repeatn=100,family="glmm_poisson",check_spa_cor=F,use_pforeach=TRUE,perm=TRUE,nrep=10000)
{
    
    if(flag_envs).ee.append("cont_selec",environment());
    
    res_contract=contract(data00,edge_explain=edge_explain,edge_cor=edge_cor,target=target);        
    group=res_contract$group;
    data1=res_contract$data1;
    data0=res_contract$data0;
    target=res_contract$target;

    data1=scale_range(data1,minvalue=0.0,maxvalue=1.0);

    
    res=selec_bestsub_stepwise(data0,data1,group,pos_x,pos_y,
                edge_param_number=edge_param_number,
                repeatn=repeatn,family=family,check_spa_cor=check_spa_cor,target=target, use_pforeach=use_pforeach,perm=perm,nrep=nrep);
    return(res);
}

#' @export
scale_range <- function(data1,minvalue=0.0,maxvalue=1.0){
    for(i in 2:ncol(data1)){
        data1[,i]=((maxvalue-minvalue)*(data1[,i]-min(data1[,i]))/(max(data1[,i])-min(data1[,i])))+minvalue;
    }
    return(data1);
}
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/ FUNCTIONS FOR CONTRACTION OF VARIABLES _/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#

##' This function converts data0 into data1, by 1) grouping variables of high correlation, and by 2) replacing each group's variables with its PCA scores. The column specified by "target" is treated as the explained variable, named "y", which is untatched.
##' @title contracts explanatory variables of high absolute correlation
##' @param data0 data.frame : explained- and explanatory-variables.
##' @param target integer or character : the column specified by "target" is treated as the explained variable "y", which is untatched. 
##' @param edge_explain real number : threshold for the variance ratio to be explained. The number of principal components "np" for representing the contracted group is determined so that [variance by the np principal components]/[total variance] exceeds "edge_explain". 
##' @param edge_cor real number : threshold for clustering; variables of correlation no less than this threshold are clustered together.
##' @param np integer : number of principal components to represent the contracted group, if assigned.
##' @return list(data1,group,data0)), where "data1"(data.frame) is contracted data. "group" (list(gid,gid_data1,ngid,vname1_orig,vname0_orig)) contains information about contraction of "data0" into "data1". "data0" is the data before contraction.
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' con=contract(data,edge_cor=0.9,edge_explain=0.6,target="Horsepower")
##'
##' con$target;
##' print(data.frame(uncontracted.name=colnames(con$data0)));
##' print(data.frame(contracted.name=colnames(con$data1),members=con$group$vname1_orig));
##' @export
contract <- function(data0,edge_explain=-1,edge_cor,np=-1,target=1){
    if(flag_envs).ee.append("contract",environment());

    if(target!=1){
        if(is.numeric(target)){
            target_id=target;
        }
        if(is.character(target)){
           target_id=which(colnames(data0)==target);            
        }

        if(target_id>1)data0=data0[,c(target_id,1:(target_id-1),(target_id+1):ncol(data0))];
        
    }
    target=colnames(data0)[1];

    colnames(data0)[1]="y";
    
    gid= assign_gid(data0,edge_cor=edge_cor);
    gid_data1=rep(0,1);    
    data1=data0;
    vname1_orig=c(colnames(data0)[1]);
        
    count_contracted_group=0;
    count_var1=2;

    for(i in 1:max(gid)){
        list0=which(gid==i);
        data0s=data0[,list0];
        
        if(length(list0)>1){
            vname1_buf=paste(colnames(data0)[list0],collapse=" + ");
            
            cat("\n   gid:",i," contracted group:  ", vname1_buf,"\n");
            
            data0s=normalize(data0s);
            if(np<0)pe=get_pcom(data0s,edge_explain=edge_explain);
            if(np>0){
                if(length(list0)==2)pe=get_pcom(data0s,edge_explain=edge_explain,np=1);
                if(length(list0)>2)pe=get_pcom(data0s,edge_explain=edge_explain,np=np);
            }
            npp=pe[[1]];
            data0s_pcomp=pe[[2]];
            if(length(list0)>2)count_contracted_group=count_contracted_group+1;

            
            for(k in 1:npp){
                if(npp>1)data1[,count_var1]=data0s_pcomp[,k];
                if(npp==1)data1[,count_var1]=data0s_pcomp;
                colnames(data1)[count_var1]=vname1_buf;
                vname1_orig=append(vname1_orig,vname1_buf);
                
                if(length(list0)>2){
                    colnames(data1)[count_var1]=sprintf("Cont.var%d.%d",count_contracted_group,k);
                    
                }
                cat((count_var1-1)," gid:",i,"  ",colnames(data1)[count_var1],"\n");
                gid_data1=append(gid_data1,i);
                count_var1=count_var1+1;
            }
            
            
        }else{
            
            data1[,count_var1]=data0s;
            colnames(data1)[count_var1]=colnames(data0)[gid==i];
            vname1_orig=append(vname1_orig,colnames(data1)[count_var1]);
            cat((count_var1-1)," gid:",i,"  ",colnames(data1)[count_var1],"\n");
            gid_data1=append(gid_data1,i);
            count_var1=count_var1+1;
        }
    }
    
    ngid=rep(0,length(gid_data1));
    for(i in 1:length(ngid))ngid[i]=sum(gid==gid_data1[i]);

    data1=data1[,1:(count_var1-1)];
    colnames(data1)=swap2dot(colnames(data1));
    vname0_orig=colnames(data0);
    colnames(data0)=swap2dot(colnames(data0));
    group=list(gid=gid,gid_data1=gid_data1,ngid=ngid,vname1_orig=vname1_orig,vname0_orig=vname0_orig);
    return(list(data1=data1,group=group,data0=data0,target=target));
}

swap2dot <- function(x){
x=gsub(",","..",x);
x=gsub("\\(",".",x);
x=gsub(")",".",x);
x=gsub("-",".",x);
x=gsub(" ",".",x);
x=gsub("\\+",".",x);
return(x);
}


agid_ref <- setRefClass("agid", fields = list(count="numeric",gid = "numeric",flag = "numeric"));

##' This function clusters explanatory variables of high absolute correlation and assign group ids for them. 
##' @title clusters explanatory variables of high absolute correlation
##' @param data0 data.frame : explained variable (1st column), and explanatory variables (2nd-last columns) to be contracted
##' @param edge_cor real number : threshold for grouping; variables of correlation no less than this threshold are grouped together
##' @return gid :  array of integer :  group-ids for explanatory variables in "data0" 
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' assign_gid(data,edge_cor=0.9);
##' 
#' @export
assign_gid <- function(data0,edge_cor){
    ndata=ncol(data0);
    agid=agid_ref$new(count=0,flag=rep(0,ndata),gid=rep(0,ndata));    
  
    for(i in 2:ndata){
        if(agid$flag[i]==0){
            agid$count = agid$count + 1;
            agid$gid[i] = agid$count;
            agid$flag[i]=1;
            waver(i,agid,data0,edge_cor)           
        }
    }

  
    return(agid$gid);
}


waver<-function(i,agid,data0,edge_cor){
   for(j in 2:ncol(data0)){
     if((i!=j)&&(agid$flag[j]==0)){
       dis=cor(data0[,i],data0[,j]);

      if(abs(dis)>edge_cor){        
        agid$flag[j] = 1;
        agid$gid[j] =  agid$count;
        waver(j,agid,data0,edge_cor);
      }
    }
  }
}

##' This function returns PCA scores for the 1th,...,"np"th principal components of "data0s" so that [explained variance]/[total variance] is larger than "edge_explain" 
##' @title gets PCA scores from data subset
##' @param data0s : data.frame or matrix
##' @param edge_explain  real number : "np" is chosen so that [explained variance]/[total variance] > "edge_explain", when "np" is not given
##' @param np  integer : number of principal components extracted
##' @param scale  TRUE/FALSE : if TRUE, variables are scaled in advance so that mean=0 and sd=1. 
##' @return list(np, data0s_pcom), where "np" (integer) is number of principal components, and "data0s_pcom" (matrix) contains PCA scores for the corresponding principal components.
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' gid=assign_gid(data,edge_cor=0.8);
##' print(gid);
##' pcom=get_pcom(data[,gid==2],edge_explain=0.85,scale=TRUE);
#' @export
get_pcom<-function(data0s,edge_explain,np=-1,scale=FALSE){
    if(flag_envs).ee.append("get_pcom",environment());
    
    pppp=(prcomp((data0s),scale=scale));
   
    varc=cumsum((pppp$sdev)^2)/sum((pppp$sdev)^2)
    if(np<0)np=min(which(varc>edge_explain));
    
      cat("    Number of principal components used: np=",np,"(",100*varc[np],"% is explained )\n")
   
   return(list(np=np,pcom=pppp$x[,1:np]));
}

normalize<-function(data){
    for(j in 1:ncol(data)){        
        data[,j]=(data[,j]-mean(data[,j]))/sd(data[,j]);
    }
    return(data);
}

      
##' This function varies "edge_cor" and counts the number of explanatry variables in "data1" that is contracted data of "data0".
##' 
##' @title counts numbers of variables in data after contraction under different values for "edge_cor"
##' @param data0 : data.frame or matrix for data.
##' @param start_edge_cor  real number : initial value for "edge_cor".
##' @param end_edge_cor  real number : last value for "edge_cor".
##' @param ndiv  integer : number of examined values for "edge_cor" between "start_edge_cor" and "end_edge_cor".
##' @param edge_explain  real number : see arguments for function "get_pcom".
##' @param show  TRUE/FALSE :  if TRUE, result is plotted.
##' @return data.frame(edge_cor,edge_explain,nvar,nvar_single)
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' check_contract_edge(data,start_edge_cor=0.2,end_edge_cor=0.95,ndiv=10,edge_explain=0.8,show=1)
##' 
#' @export
check_contract_edge <- function(data0,start_edge_cor=0.2,end_edge_cor=0.7,ndiv=10,edge_explain,show=TRUE){
    result=data.frame(matrix(0,ncol=4,nrow=ndiv));
    colnames(result)=c("edge_cor","edge_explain","nvar","nvar_single");
    
    for(i in 1:ndiv){
        edge_cor=start_edge_cor+i*(end_edge_cor-start_edge_cor)/as.double(ndiv);
        res_con=contract(data0,edge_cor=edge_cor,edge_explain=edge_explain);
        
        nvar=length(res_con$group$gid_data1)-1;
        ngid=res_con$group$ngid;
        nvar_single=sum(ngid[2:length(ngid)]==1);
        result[i,]=c(edge_cor,edge_explain,nvar,nvar_single);
    }
    if(show){
        plot(result$edge_cor,result$nvar,col="blue",type="b",xlab="edge_cor",ylab="number of explanatory variables",ylim=c(0,max(max(result$nvar),max(result$nvar_single))));
        points(result$edge_cor,result$nvar_single,col="red",type="b")
        legend("topleft", legend = c("total","uncontracted"), col = c("blue","red"),cex=0.9,pch=1);
        }
    return(result);
}




#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/ FUNCTIONS FOR BESTSUBSET-STEPWISE MODEL SELECTION  _/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#


##' This function conducts best subset model selection, permutation tests and stepwise model selection. Whether explanatory variables in the best model have statistically strong effect on the explained variable is examined. 
##'
##' Details.
##' @title conducts best subset model selection, permutation test and stepwise model selection.
##' @param data0 data.frame : uncontracted data
##' @param data1 data.frame : constracted data
##' @param group (list(gid,gid_data1,ngid,vname1_orig,vname0_orig)) : information about contraction of "data0" into "data1".
##' @param pos_x array of numeric : x-coordinates of spatial position for each sampled data in "data00".
##' @param pos_y array of numeric : y-coordinates of spatial position for each sampled data in "data00".
##' @param edge_param_number real number: this parameter constrains the number of explanatory variables in the subset model selection and stepwise model selection, so that [sample size]/[number of free parameter] >= "edge_param_number".
##' @param repeatn integer : the number of model evaluations proccessed as a single job for each CPU core (This parameter is meaningful only when "use_pforeach" is TRUE).
##' @param family character : error distribution. Possible choices are "glmm_poisson", "poisson","gaussian".
##' @param check_spa_cor TRUE/FALSE : if TRUE, the best model with residuals significantly correlated with c(pos_x, pos_y) are removed repeatedly until the correlation becomes non-significant.
##' @param use_pforeach  TRUE/FALSE : if TRUE, pforeach is used instead of foreach in best subset selection and permutation tests (when "perm" is TRUE).
##' @param perm  TRUE/FALSE : if TRUE, permutation tests for explained variables are conducted.
##' @param nrep  integer : number of resampling in the permutation test.
##' @param target character : arbitrary name for explained variable "y".
##' @return list(target=target,result=result,data0=data0,data1=data1,group=group,
##' family=family,bestmodels=res_sb$bestmodels,bestmodel_stepwise=res_ss$bestmodel,
##' bestmodel=res_sb$bestmodel,pos_x=pos_x,pos_y=pos_y);
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' con=contract(data,edge_cor=0.9,edge_explain=0.6,target="Horsepower");
##' data1=con$data1;
##' for(i in 2:ncol(data1))data1[,i]=(data1[,i]-min(data1[,i]))/(max(data1[,i])-min(data1[,i]));
##'
##' res=selec_bestsub_stepwise(con$data0,data1,con$group,
##' edge_param_number=3,repeatn=50,family="gaussian",
##' target=con$target, use_pforeach=FALSE,perm=TRUE);
##'
##' print(res);
##' summary(res);
#' @export
selec_bestsub_stepwise <- function(data0,data1,group,pos_x=NULL,pos_y=NULL,
                    edge_param_number=3,
                    repeatn=100,family="glmm_poisson",check_spa_cor=FALSE,target="",use_pforeach=TRUE,perm=TRUE,nrep=10000)
{

    if(flag_envs).ee.append("selec_bestsub_stepwise",environment());

    
    print("best subset selection...");
    if(use_pforeach)cat("using pforeach\n");
    if(!use_pforeach)cat("using foreach\n");

    cat("[sample size]/[free parameters] <= ",edge_param_number,"\n");

  
    res_sb=selec_bestsub(data0,data1,group,pos_x=pos_x,pos_y=pos_y,edge_param_number=edge_param_number,repeatn=repeatn,family=family,check_spa_cor=check_spa_cor,target=target,use_pforeach=use_pforeach,perm=perm,nrep=nrep);

    result=res_sb$result;
    
        
    print("stepwise selection...");
    res_ss=selec_stepwise(result,data0,edge_param_number,family=res_sb$bestmodel$family,scale01=TRUE);
   
    result=res_ss$result;
    result=put_sigmark(result);
    cat("Best model by best-subset-selection:\n");
    print(res_sb$bestmodel$form_name_coef);
    cat("Best model by stepwise-selection:\n");
    print(res_ss$bestmodel$form_name_coef);
    cat("statistically strong factors:\n\n");

 
    print_result(result,res_sb$bestmodel$family);

    print("number of models:");
    print(res_sb$n_model_total);
    
    print("number of models fitted:");
    print(res_sb$n_model);

    if(check_spa_cor){
        print("number of removed best models due to spatial correlation:");
        print(length(res_sb$bestmodels_spa_cor$list_remove));
    } 


    outer_bests=cbind(sigmark=result$sigmark,res_sb$outer_bests);
    
    result_all=list(target=target,result=result,data0=data0,data1=data1,group=group,family=family,bestmodels=res_sb$bestmodels,bestmodel_stepwise=res_ss$bestmodel,bestmodel=res_sb$bestmodel,pos_x=pos_x,pos_y=pos_y,bestmodels_spa_cor=res_sb$bestmodels_spa_cor,outer_bests=outer_bests);

    return(structure(result_all,class="selec_bestsub_stepwise"));
}


    
##' This function shows summary of the result of "cont_selec" or "selec_bestsub_stepwise". See "cont_selec" for an example.
##' @title shows summary of "cont_selec" or "selec_bestsub_stepwise" result.
##' @param res : output of "selec_bestsub_stepwise"
##' @author Hiroshi C. Ito
#' @export
print.selec_bestsub_stepwise <- function(res){
    print_result(res$result,res$family);
}

##' This function shows summary of the result of "cont_selec" or "selec_bestsub_stepwise". See "cont_selec" for an example.
##' @title shows summary of "cont_selec" or "selec_bestsub_stepwise" result.
##' @param res : output of "selec_bestsub_stepwise"
##' @author Hiroshi C. Ito
#' @export
summary.selec_bestsub_stepwise <- function(res){
    print_result(res$result,res$family);
}


##' This function puts sigmarks "*","**" and "***" for p-values <0.05, <0.01 and <0.001 to the result output of "selec_stepwise".
##' @title puts sigmarks
##' @param result : output of "selec_bestsub_stepwise"
##' @author Hiroshi C. Ito
#' @export
put_sigmark <- function(result){
    for(i in 1:length(result$sigmark)){
        result$sigmark[i]="   ";
        if(result$nmodel[i]*result$stepwise[i]==result$nmodel[1]){
            if(result$pv[i]<0.05)result$sigmark[i]="  *"
            if(result$pv[i]<0.01)result$sigmark[i]=" **"
            if(result$pv[i]<0.001)result$sigmark[i]="***"
        }
    }
    return(result);
}



print_result <- function(result,family){
    show_ratio=TRUE;
    if(family=="gaussian")show_ratio=FALSE;
    for(i in 2:nrow(result)){
        if(abs(result$nmodel[i])>0){
            cat(result$sigmark[i],result$contname[i],":",result$vname0_maxcor[i],"\n     ",result$nmodel[i],"/",result$nmodel[1], "  step:",result$stepwise[i],"  pvalue:",result$pv[i]," coef:",result$coef[i], "(se: ",result$coef_se[i],")");
            if(show_ratio==TRUE)cat(" ratio:",exp(result$coef[i]));
            cat("\n");
        }
    }
}


#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/ FUNCTIONS FOR BEST-SUBSET MODEL SELECTION   _/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#



##' This function conducts best-subset model selection
##'
##' Details.
##' @title best-subset model selection
##' @param data0 data.frame : uncontracted data
##' @param data1 data.frame : constracted data
##' @param group (list(gid,gid_data1,ngid,vname1_orig,vname0_orig)) : information about contraction of "data0" into "data1".
##'
##' @param pos_x array of numeric : x-coordinates of spatial position for each sampled data in "data0".
##' @param pos_y array of numeric : y-coordinates of spatial position for each sampled data in "data0".
##' @param edge_param_number real number: this parameter constrains the number of explanatory variables in the subset model selection and stepwise model selection, so that [sample size]/[number of free parameter] >= "edge_param_number".
##' @param repeatn integer : the number of model evaluations proccessed as a single job for each CPU core (This parameter is meaningfur only when parallel proccessing is conducted by library "pforeach").
##' @param family character : error distribution. Possible choices are "glmm_poisson", "poisson","gaussian".
##' @param check_spa_cor TRUE/FALSE : if TRUE, the best model with residuals significantly correlated with c(pos_x, pos_y) are removed repeatedly until the correlation becomes non-significant. 
##' @param target character : arbitrary name for explained variable "y".
##' @param use_pforeach  TRUE/FALSE : if TRUE, pforeach is used instead of foreach.
##' @param perm  TRUE/FALSE : if TRUE, permutation tests for explained variables are conducted.
##' @param nrep  integer : number of resampling in the permutation test.
##' @param daic  numeric : models with Delta-AIC no more than "daic" are examined. Default value is 2.0.
##' @return list(result=result,bestmodels=aa,bestmodel=bestmodel,
##' n_model=n_model,n_model_total=n_model_total)
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' con=contract(data,edge_cor=0.9,edge_explain=0.6,target="Horsepower");
##' data1=con$data1;
##' for(i in 2:ncol(data1))data1[,i]=(data1[,i]-min(data1[,i]))/(max(data1[,i])-min(data1[,i]));
##'
##' res=selec_bestsub(con$data0,data1,con$group,edge_param_number=3,repeatn=50,family="gaussian",target=con$target, use_pforeach=FALSE,perm=TRUE,daic=2.0);
##'
##' print(res$result);
##'
##' 
##' 
#' @export
selec_bestsub <- function(data0,data1,group,pos_x=NULL,pos_y=NULL,
                    edge_param_number=3,
                    repeatn=100,family,check_spa_cor=FALSE,target="",use_pforeach=TRUE,perm=TRUE,daic=2.0,check_spa_cor_nmodel=50,nrep=10000){

    if(flag_envs).ee.append("selec_bestsub",environment());
    
    res_selec_bestsubset=selec_bestsubset(data1,daic=-1,repeatn=repeatn,edge_param_number=edge_param_number,family=family,use_pforeach=use_pforeach);

    aa0=res_selec_bestsubset$dataframe;

    check_spa_cor_nmodel=min(check_spa_cor_nmodel,nrow(aa0));
    spa_cor=NULL;
    if(check_spa_cor){
        print(paste("checking spatial correlation in residuals for best", check_spa_cor_nmodel, "models..."));
        spa_cor=get_spa_cor(aa0[1:check_spa_cor_nmodel,],data1,pos_x,pos_y);        
        if(length(spa_cor$list_remove)){
            ##aa0=aa[-spa_cor$list_remove,];
            aa0=aa0[-spa_cor$list_remove,];
            aa0$delta_aic=aa0$aic-aa0$aic[1];
        }
    }
    
    nbest=sum(aa0$delta_aic<=daic);
    aa=aa0[1:nbest,];
    
    
    print_bestsubset(aa);
    
    bestmodel=get_best_subset(aa[1,],data1);
    
    result=make_result(bestmodel,aa,data0,data1,group,use_pforeach=use_pforeach,perm=perm,nrep=nrep);


    outer_bests=result[,2:6];
    outer_bests[,2:5]=outer_bests[,2:5]*0;
    colnames(outer_bests)=c("contname","delta_aic","aic","order","form_name");
    for(i in 2:nrow(outer_bests)){
        if(result$pv[i]>-0.5){
            myname=result$contname[i];
            if(sign(result$coef[i])>0)myname=paste0("+",myname);
            if(sign(result$coef[i])<0)myname=paste0("-",myname);
            
            have=grep(myname,aa0$form_name);
            ii=which(have-c(1:length(have))>0)[1];
            outer_bests$delta_aic[i]=aa0$delta_aic[ii];
            outer_bests$aic[i]=aa0$aic[ii];
            outer_bests$form_name[i]=aa0$form_name[ii];
            outer_bests$order[i]=ii;
            
            }

        }
    
    
    n_model=length(aa0$aic);
    n_model_total=res_selec_bestsubset$total;

    
    return(list(result=result,bestmodels=aa,bestmodel=bestmodel,n_model=n_model,n_model_total=n_model_total,bestmodels_spa_cor=spa_cor,outer_bests=outer_bests));
}




selec_bestsubset <- function(data,daic=-1,repeatn=50,edge_param_number=3, family,use_pforeach=TRUE){
    if(flag_envs).ee.append("selec_bestsubset",environment());    
 
    
    data = subset(data, complete.cases(data));    
    nvv = ncol(data)-1;                     
    n = (2^nvv) -1;
    if(nvv>32)cat("Too many explanatory variables!!!\n");
    vcomb=(sapply(1:n,intToBits)[1:nvv,]>0);
    vcomb=t(vcomb);
   
    depname = colnames(data)[1];
    vname = colnames(data)[2:ncol(data)];
    ##data_buf=transform(data,id=c(1:nrow(data)));
    data_buf=cbind(data,id=seq(nrow(data)));
    
    nsample=nrow(data);        
    nparam=rowSums(vcomb);
    plimit=as.integer(nsample/edge_param_number)-1;
    lists=which(nparam<=plimit);
    lists_glmm=which(nparam<=plimit-1);
    vcomb1=vcomb[lists,];
    n1=length(lists);
    vcomb1_glmm=vcomb[lists_glmm,];
    n1_glmm=length(lists_glmm);
    
    reg_param=list(n=n1,data=data_buf,vname=vname,vcomb=vcomb1,depname=depname,family=family,repeatn=repeatn,edge_param_number=edge_param_number);
    if(family=="glmm_poisson")reg_param$family="glmm_poisson_fixs";
    
    reg_param_glmm=list(n=n1_glmm,data=data_buf,vname=vname,vcomb=vcomb1_glmm,depname=depname,family="glmm_poisson",repeatn=repeatn,edge_param_number=edge_param_number);

    cat("\n fitting",n1,"models...\n");
    ni= 1+as.integer((n1-1)/repeatn);
    if(use_pforeach){
        results=pforeach::pforeach(i=1:ni, .combine='rbind',.errorhandling="stop")({    
            regressions(i,reg_param);     
        });
    }else{
        results=foreach::foreach(i=1:ni,.combine='rbind',.packages="foreach")%do%{    
            regressions(i,reg_param);     
        };
    }
    

    if(family=="glmm_poisson"){
        cat("\n fitting",n1_glmm,"models with random effect...\n");
        ni= 1+as.integer((n1_glmm-1)/repeatn);
               
        if(use_pforeach){
            results_glmm=pforeach::pforeach(i=1:ni, .combine='rbind',.errorhandling="stop")({    
                regressions(i,reg_param_glmm);     
            });
        }else{
            results_glmm=foreach::foreach(i=1:ni,.combine='rbind',.packages="foreach")%do%{    
                regressions(i,reg_param_glmm);     
            };
        }
        
        
        form_name=append(results$form_name,results_glmm$form_name);
        aic=as.numeric(append(results$aic,results_glmm$aic));
        form=(append(results$form,results_glmm$form));
    }else{
        
        form_name=results$form_name;
        aic=as.numeric(results$aic);
        form=results$form;
    }
    

    if(family=="glmm_poisson"){
        total=2*n;
    }else{
        total=n;
    }
    
    oo<-order(aic);
    aic=aic[oo];
    delta_aic=aic-aic[1];
    form=form[oo];
    form_name=form_name[oo];
    
    aa0=data.frame(aic=aic, delta_aic=delta_aic, form=form,form_name=form_name,stringsAsFactors=FALSE);
    if(daic<0){
        return(list(dataframe=aa0,total=total));
    }else{
        nout=which(delta_aic<=daic);
        return(list(dataframe=aa0[nout,],total=total));
    }
}


print_bestsubset <- function(result,model_number=20)
{
  cat("\n Delta-AIC       AIC         Formula\n");  
  for(i in 1:nrow(result)){                       
      cat(sprintf("%.5f   %.5f  %s\n", result$delta_aic[i], result$aic[i], result$form_name[i]));
  }
}


#' @export
regressions<-function(k,reg_param){
    n=reg_param$n;
    repeatn=reg_param$repeatn;
    edge_param_number=reg_param$edge_param_number;

    nj=repeatn*k;
    if(nj>n)nj=n;

    if(k%%40==0)cat(sprintf("%d%%",as.integer(100.0*repeatn*(k-1)/n)));
    
    res=foreach::foreach(j=((k-1)*repeatn+1):nj,.combine='rbind',.packages="foreach")%do%{    
        regression(j,reg_param);   
    };
    
    return(res);    
}

aic_init=100000;


regression <- function(i,reg_param){
   ## if(flag_envs).ee.append("regression",environment());
    n=reg_param$n;
    data=reg_param$data;
    vname=reg_param$vname;
    depname=reg_param$depname;
    vcomb=reg_param$vcomb;
    family=reg_param$family;
    edge_param_number=reg_param$edge_param_number;
    nsample=nrow(data);
        
    nparam=sum(vcomb[i,]);
    
    form = reformulate(vname[vcomb[i,]], depname);
    

    aic=aic_init;
    form_name="";
      
    formc=form2char(form);

 ##   plimit=as.integer(nsample/edge_param_number)-1;
    
    if(family=="glmm_poisson"){
        ##if(nparam<=(plimit-1)){
      ## same with    if(nsample>=(nparam+1)*edge_param_number)
            gg_poi = glmmML::glmmML(form,data=data,family=poisson,cluster=id,start.sigma=0.001);
            aic=gg_poi$aic;
            form_name=get_form_name(gg_poi$coefficients,family);            
        ##}
    }
    if(family=="glmm_poisson_fixs"){
        ##if(nparam<=plimit){
            gg_poi = glmmML::glmmML(form,data=data,family=poisson,cluster=id,fix.sigma=TRUE,start.sigma=0.0);
            aic=gg_poi$aic;
            form_name=get_form_name(gg_poi$coefficients,"glmm_poisson_fixs");
        ##}
    }
    
    if((family=="gaussian")||(family=="poisson")){
       ## if(nparam<=plimit){
            gg_poi=glm(form,data=data,family=family);
            aic=gg_poi$aic;
            
            form_name=get_form_name(gg_poi$coefficients,family);
        ##}
    }
    
         
    res_unit=data.frame(form_name=form_name,aic=aic,form=formc,stringsAsFactors=FALSE);

    return(res_unit);	
}


form2char <- function(form){
        temp = as.character(form);    
        formc=paste(temp[2], "~", temp[3]);
        return(formc);
}

get_form_name <- function(coef,family){
    coef=coef[2:length(coef)];
    cnames=names(coef);
    for (k in 1:length(cnames)){
        if(coef[k]>=0)cnames[k]=paste("+",cnames[k],sep="");
        if(coef[k]<0)cnames[k]=paste("-",cnames[k],sep="");
    }
    form_name=paste(cnames,collapse=" ");
    form_name=paste("y(",family,") ~ ",form_name,sep="");
    return(form_name);
}



get_family <- function(form_name){
    buf=(as.character(as.formula(form_name)))[2];    
    buf=gsub(")","",buf);
    buf=gsub("y\\(","",buf);
    return(buf);
}

get_best_subset <- function(aa_top,data1){

    if(flag_envs).ee.append("get_best_subset",environment());
    
    ##data1_buf=transform(data1,id=c(1:nrow(data1)));
    data1_buf=cbind(data1,id=seq(nrow(data1)));
    
    family=get_family(aa_top$form_name);
    family_name=family;
    if((family=="glmm_poisson")||(family=="glmm_poisson_fixs")){
        if(family=="glmm_poisson"){
            fix.sigma=FALSE;start.sigma=0.001;
        }
        if(family=="glmm_poisson_fixs"){
            fix.sigma=TRUE;start.sigma=0.0;
            family_name="glmm_poisson_fixs";
        }
        
        g3=glmmML::glmmML(as.formula(aa_top$form),family=poisson,data=data1_buf,cluster=id,fix.sigma=fix.sigma,start.sigma=start.sigma);
                
        coef = g3$coefficients
        se = g3$coef.sd
        pv0=signif(1 - pchisq((coef/se)^2, 1));
      
        print(glmmML::summary.glmmML(g3));        
    }else{
        if(family=="gaussian"){
            g3=glm(as.formula(aa_top$form),family="gaussian",data=data1);            
        }
        if(family=="poisson"){
            g3=glm(as.formula(aa_top$form),family="poisson",data=data1);            
        }
        print(summary(g3));
        pv0=(summary(g3)$coefficients)[,4];
        se=(summary(g3)$coefficients)[,2];
    }
       
    coef=g3$coefficients;
    nmodel=array(0,length(coef));
    
    form_name_coef=get_form_name_coef(g3$coefficients,family);
    
    resi=get_resi(g3,data1,family=family);
    r2=1.0-var(resi)/var(data1[,1]);
    
    return(list(pv_coef=data.frame(pv=pv0,coef=coef,coef_se=se,nmodel=nmodel),model=g3,form_name_coef=form_name_coef,r2=r2,family=family_name));
}


get_resi <- function(gg,data1,family){  
    coef=gg$coefficients;
    return(get_resi0(coef,data1,family=family));

}

get_resi0  <- function(coef,data1,family){
    yy=data1[,1];
    y_est=yy*0;
    y_est=y_est+coef[1];
    vname=names(coef);
    for(i in 2:length(coef)){
        value=data1[,which(colnames(data1)==vname[i])]
        y_est=y_est+value*coef[i]; 
    }
    if((family=="glmm_poisson")||(family=="glmm_poisson_fixs")||(family=="poisson")){
        resi=yy-exp(y_est);
    }
    if(family=="gaussian"){
        resi=yy - y_est;
        }
    return(as.vector(resi));
}



#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/ FUNCTIONS FOR CHECKING SPATIAL CORRELATION _/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#

get_spa_cor <- function(aa,data,pos_x,pos_y){

    naa=length(aa$aic);
    list_spa_cor=c(NULL);
    spa_cor=data.frame(moran=rep(1.0,naa),gear=rep(1.0,naa));
    ##data1_buf=transform(data,id=c(1:nrow(data)))
    data1_buf=cbind(data,id=seq(nrow(data)));
    for(i in 1:naa){
        family=get_family(aa$form_name[i]);
        
        if(family=="glmm_poisson")gg=glmmML::glmmML(as.formula(aa$form[i]),family=poisson,data=data1_buf,cluster=id,start.sigma=0.001);
        if(family=="glmm_poisson_fixs")gg=glmmML::glmmML(as.formula(aa$form[i]),family=poisson,data=data1_buf,cluster=id,fix.sigma=TRUE,start.sigma=0.0);
        
        if((family=="gaussian")||(family=="poisson")){
            gg=glm(as.formula(aa$form[i]),family=family,data=data1_buf);
        }
        
        resi=get_resi(gg,data,family=family);
        spac=check_spatial_cor(resi,pos_x,pos_y);
        spa_cor$moran[i]=spac[1];
        spa_cor$gear[i]=spac[2];
        if(min(spac)<0.05)list_spa_cor=append(list_spa_cor,i);
        
    }

    list_remove=integer(0);
    if(length(list_spa_cor)>0){
    list_remove=list_spa_cor*0;
    for(i in 1:length(list_spa_cor)){
        if(list_spa_cor[i]==i)list_remove[i]=list_spa_cor[i];
    }
    list_remove=list_remove[which(list_remove>0)];
    }
    return(list(list_remove=list_remove,list_sig=list_spa_cor,pvalue=spa_cor));
}
  

check_spatial_cor <- function(resi,pos_x,pos_y,show=0){
     
    nsamp=length(resi);
    pos=as.matrix(cbind(pos_x,pos_y));
    postri=spdep::tri2nb(pos,row.names=as.character(1:nsamp));
   
    mm=spdep::moran.test(resi, spdep::nb2listw(postri, style ="W"),alternative="greater");
    if(show==1)print(mm);
    mg=spdep::geary.test(resi, spdep::nb2listw(postri, style ="W"),alternative="greater");
    if(show==1)print(mg);
        
    return(c(mg$p.value,mm$p.value));
}


#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/  FUNCTIONS FOR POST-PROCESSING FOR selec_bestsubset _/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#

make_result <- function(bestmodel,aa,data0,data1,group,perm=TRUE,use_pforeach=TRUE,nrep=nrep){
    if(flag_envs).ee.append("make_result",environment());
    pv_coef=bestmodel$pv_coef;
    pv_param=pv_coef$pv;
    if(perm){
        print("permutation test...");
        pv_coef_perm=get_perm_pv_coef(aa[1,],pv_coef,data1,use_pforeach=use_pforeach,nrep=nrep);
    }else{
        pv_coef_perm=pv_coef;
    }
    pv_coef_perm_nmodel=count_nmodel(pv_coef_perm,aa);

    
    pv_coef1=adjust_pv_coef(pv_coef_perm_nmodel,colnames(data1),pv_param);

    ncol1=ncol(data1);
    result=data.frame(sigmark=rep("",ncol1),
                      contname=colnames(data1),
                      nmodel=pv_coef1$nmodel,
                      stepwise=rep(0,ncol1),
                      pv=pv_coef1$pv,
                      pv_paramet=pv_coef1$pv_paramet,
                      coef=pv_coef1$coef,
                      coef_se=pv_coef1$coef_se,
                      vname0_maxcor=rep("",ncol1),
                      cor=rep(0.0,ncol1),
                      gid=group$gid_data1,
                      ngid=group$ngid,
                      stringsAsFactors=FALSE);
    rownames(result)=as.character(1:nrow(result));

    result=add_maxcor_vname0(result,data1,group,data0);
    
return(result);
}

get_perm_pv_coef <- function(aa_top,pv_coef,data1,use_pforeach=TRUE,nrep){
    ##if(flag_envs).ee.append("get_perm_pv_coef",environment());
    ni=length(pv_coef$pv);
    if(use_pforeach){
        pv_perm=pforeach::pforeach(i=2:ni,.combine='c',.errorhandling="stop")({
            perm_each(i,pv_coef,aa_top,data1,nrep=nrep);
        });
    }else{
        pv_perm=foreach::foreach(i=2:ni,.combine='c',.packages="foreach")%do%{
            perm_each(i,pv_coef,aa_top,data1,nrep=nrep);
        }        
    }

     pv_coef$pv[2:ni]=pv_perm;

  return(pv_coef);
}


#' @export
perm_each <- function(i,pv_coef,aa_top,data1,seed=52,nrep=10000){
    pv0=pv_coef$pv;
    coef=pv_coef$coef;
    names(pv0)=rownames(pv_coef);
    form=as.formula(aa_top$form);
    family=get_family(aa_top$form_name);
    if(family=="glmm_poisson_fixs")family1="poisson";
    if(family=="glmm_poisson")family1="quasipoisson";
    if(family=="poisson")family1="poisson";
    if(family=="gaussian")family1="gaussian";
    
    perm=glmperm::prr.test(form,var=names(pv0)[i],data=data1,family=family1,nrep=nrep,seed=seed);
    ppv=perm$p.value.perm$p0[[1]];
    print(paste("p-value of",names(pv0)[i]," : ",ppv,",  se:",perm$p.value.perm.se$se.p0));
   
  return(ppv);
}



count_nmodel  <- function(pv_coef,aa){
    if(flag_envs).ee.append("count_nmodel",environment());
    
    taf=aa$form_name;
    pv_coef_buf=pv_coef;
    pv_coef_buf$nmodel[1]=length(taf);
    coef=pv_coef$coef;
    for(i in 2:nrow(pv_coef)){
        ta=rownames(pv_coef)[i];
        if(coef[i]>0)ta=paste("\\+",ta,sep="");
        if(coef[i]<=0)ta=paste("-",ta,sep="");
        
        ta_count=count_target_var(ta,taf);
        if(coef[i]>0)pv_coef_buf$nmodel[i]=ta_count;
        if(coef[i]<0)pv_coef_buf$nmodel[i]=-1*ta_count;
    }
    return(pv_coef_buf);
}

count_target_var <- function(ta,taf){
    return(length(grep(ta,taf)));      
}


adjust_pv_coef <- function(pv_coef,vnames,pv_paramet){
    if(flag_envs).ee.append("adjust_pv_coef",environment());
    pv=pv_coef$pv;
    coef=pv_coef$coef;
    coef_se=pv_coef$coef_se;
    nmodel=pv_coef$nmodel;
    
    nv=length(vnames);
    coef1=rep(0.0,nv);
    coef_se1=rep(0.0,nv);
    pv1=rep(-1,nv);
    pv_paramet1=rep(-1,nv);
    
    nmodel1=rep(0,nv);
    nmodel1[1]=nmodel[1]; 
    coef1[1]=coef[1];
    
    for(j in 2:length(coef)){
        i=which(vnames==rownames(pv_coef)[j]);
        coef1[i]=coef[j];
        coef_se1[i]=coef_se[j];
        pv1[i]=pv[j];
        pv_paramet1[i]=pv_paramet[j];
        nmodel1[i]=nmodel[j];      
    }

    names(nmodel1)=vnames;
    names(pv1)=vnames;
    names(pv_paramet1)=vnames;
    names(coef1)=vnames;
    return(data.frame(pv=pv1,pv_paramet=pv_paramet1,coef=coef1,coef_se=coef_se1,nmodel=nmodel1));
}


add_maxcor_vname0 <-function(result,data1,group,data0){
        gid=group$gid;
        gid_data1=group$gid_data1;

    for(i in 2:nrow(result)){    
        cors=cor(data1[,i],data0);
        if(sum(gid==gid_data1[i])!=2){
            myname0=colnames(data0)[which.max(abs(cors))];
        }else{
            myname0=colnames(data0)[(order(-1*abs(cors)))[1:2]];
            myname0=paste(myname0,collapse=" + ");
        }
        mycor=cors[which.max(abs(cors))];
        result$vname0_maxcor[i]=myname0;
        result$cor[i]=mycor;    
    }
    return(result);
}

#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/_/_/ FUNCTIONS FOR STEPWISE SELECTION _/_/_/_/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#


##' This function conducts stepwise model selection, using the best model selected by "selec_bestsub" as the initial model. 
##'
##' @title stepwise model selection
##' @param result data.frame : result output of "selec_bestsub"
##' @param data0 data.frame : explained variable "y" and uncontracted explanatory variables.
##' @param edge_param_number real number: this parameter constrains the number of explanatory variables in the subset model selection and stepwise model selection, so that [sample size]/[number of free parameter] >= "edge_param_number".
##' @param family character : error distribution. Possible choices are "glmm_poisson", "poisson","gaussian".
##' @param scale01  TRUE/FALSE : if TRUE, explanatory variables are scaled in advance so that their ranges are from 0 to 1.
##' @return list(result=result,bestmodel=bestmodel)
##' @author Hiroshi C. Ito
##' @examples
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' con=contract(data,edge_cor=0.9,edge_explain=0.6,target="Horsepower");
##' data1=con$data1;
##' for(i in 2:ncol(data1))data1[,i]=(data1[,i]-min(data1[,i]))/(max(data1[,i])-min(data1[,i]));
##'
##' res_sb=selec_bestsub(con$data0,data1,con$group,edge_param_number=3,repeatn=50,family="gaussian",target=con$target, use_pforeach=FALSE,perm=TRUE,daic=2.0);
##' print(res_sb$result);
##' 
##' res_ss=selec_stepwise(res_sb$result,con$data0,edge_param_number=3,family=res_sb$bestmodel$family,scale01=TRUE);
##' print(put_sigmark(res_ss$result));
##'
#' @export
selec_stepwise <- function(result,data0,edge_param_number,family,scale01=FALSE){
    if(flag_envs).ee.append("selec_stepwise",environment());

    if(scale01)data0=scale_range(data0,minvalue=0.0,maxvalue=1.0);
    
    form_step0_x=paste(result$vname0_maxcor[result$pv>-0.5],collapse=" + ");
    form_step0=paste("y ~",form_step0_x);
    
    form_step_all_x=paste(colnames(data0)[2:(ncol(data0))],collapse=" + ");
    form_step_all=paste("y ~",form_step_all_x);
    
    print(form_step0);
    
    forms0=as.formula(form_step0);
    forms_all=as.formula(form_step_all);
    li=list(upper=forms_all,lower=~1);

    ##data0_buf=transform(data0,id=c(1:nrow(data0)));
    data0_buf=cbind(data0,id=seq(nrow(data0)));
    if(family=="glmm_poisson_fixs")g0=glmmML::glmmML(forms0,cluster=id,start.sigma=0.0,fix.sigma=TRUE,data=data0_buf,family=poisson);
    
    if(family=="glmm_poisson")g0=glmmML::glmmML(forms0,cluster=id,start.sigma=0.001,data=data0_buf,family=poisson);

    if((family=="gaussian")||(family=="poisson")){
        g0=glm(forms0,data=data0,family=family);
    }
    
    nsample=nrow(data0);
        
    plimit=as.integer(nsample/edge_param_number)-1;
    if(family=="glmm_poisson")plimit=plimit-1;
    
    print("Number of paremeters are limited upto:");
    print(plimit);
    g1=stepAIC_modified(g0,direction="both",scope=li,param_limit=plimit,trace=0);
        
    g1form=get_form_name(g1$coefficients,family);
    
    nmodel1_step=count_nmodel_stepwise(g1$coefficients,result);
    form_name_coef=get_form_name_coef(g1$coefficients,family);
    resi=get_resi(g1,data0,family=family);
    r2=1.0-var(resi)/var(data0[,1]);

    result$stepwise=nmodel1_step;
    bestmodel=list(model=g1,form_name=g1form,form_name_coef=form_name_coef,nmodel=nmodel1_step,r2=r2);
        
return(list(result=result,bestmodel=bestmodel));
}


count_nmodel_stepwise  <- function(coef,result){
    taf=get_form_name(coef,"");
    nmodel1_step=result$stepwise*0;
    
    names(nmodel1_step)=result$contname;
    for(i in 2:length(nmodel1_step)){
        
        if(result$ngid[i]!=2){
            ta=result$vname0_maxcor[i];
            if(length(grep(paste("\\+",ta,sep=""),taf))>0)nmodel1_step[i]=1;
            if(length(grep(paste("-",ta,sep=""),taf))>0)nmodel1_step[i]=-1;
        }else{
            ta=strsplit(result$vname0_maxcor[i],split=" ")[[1]][c(1,3)];
            nmodel1_step[i]=0;
            if(length(grep(paste("\\+",ta[1],sep=""),taf))>0)nmodel1_step[i]=nmodel1_step[i]+1;
            if(length(grep(paste("-",ta[1],sep=""),taf))>0)nmodel1_step[i]=nmodel1_step[i]-1;
            if(length(grep(paste("\\+",ta[2],sep=""),taf))>0)nmodel1_step[i]=nmodel1_step[i]+1;
            if(length(grep(paste("-",ta[2],sep=""),taf))>0)nmodel1_step[i]=nmodel1_step[i]-1;
            if(nmodel1_step[i]>0)nmodel1_step[i]=1;
            if(nmodel1_step[i]<0)nmodel1_step[i]=-1;
            
        }
    }
    
    nmodel1_step=nmodel1_step*sign(result$cor);
    return(nmodel1_step);
}

get_form_name_coef <- function(coefb,family){
    if(family=="gaussian"){
        form_name_coef=sprintf("(%s) E(y) = %.2f",family,coefb[1]);
    }else{
        form_name_coef=sprintf("(%s) E(ln(y)) = %.2f",family,coefb[1]);
    }
    for(i in 2:length(coefb)){
        plus="";
        if(coefb[i]>0)plus="+";
        
        varb=sprintf("%s%.2f[%s]",plus,coefb[i],names(coefb)[i]);
        form_name_coef=paste(form_name_coef,varb);
    }
    return(form_name_coef);
}

check_cor  <- function(data,edge_cor){
    print("correlations higher than edge_cor among explanatory variables:");
    vname=colnames(data);
    for(i in 2:(ncol(data)-1)){
        for(j in (i+1):ncol(data)){
            corv=cor(data[,i],data[,j]);
            if(abs(corv)>edge_cor){
                cat(corv,vname[i],vname[j],"\n");
                }
        }
    }
}


##' This function check interactions among statistically-strong explanatory variables in the result output of "cont_selec" 
##'
##' Details.
##' @title Checks interactions in the result output of "cont_selec"
##' @param res list : result output of "cont_selec".
##' @param edge_param_number real number: this parameter constrains the number of explanatory variables in the subset model selection and stepwise model selection, so that [sample size]/[number of free parameter] >= "edge_param_number".
##' @param repeatn integer : the number of model evaluations proccessed as a single job for each CPU core (This parameter is meaningfur only when parallel proccessing is conducted by library "pforeach").
##' @param use_pforeach  TRUE/FALSE : if TRUE, pforeach is used instead of foreach.
##' @param perm  TRUE/FALSE : if TRUE, permutation tests for explained variables are conducted.
##' @param check_spa_cor TRUE/FALSE : if TRUE, the best model with residuals significantly correlated with c(pos_x, pos_y) are removed repeatedly until the correlation becomes non-significant.
##' @param check_spa_cor_nmodel integer : for the top [check_spa_cor_nmodel] models, spatial correlation of residuals are checked. 
##' @param omit_contracted TRUE/FALSE : if TRUE, this function checks interactions among only explanatory variables that are uncontracted ones (i.e., res$group$ngid==1).
##' @param edge_cor real number: explanatory variables in result with absolute correlations higher than "edge_cor" are printed out. 
##' @return list(target=target,result=result,data0=data0,data1=data1,group=group,
##' family=family,bestmodels=res_sb$bestmodels,bestmodel_stepwise=res_ss$bestmodel,
##' bestmodel=res_sb$bestmodel,pos_x=pos_x,pos_y=pos_y);
##' @author Hiroshi C. Ito
##' @examples
##' 
##' data(Cars93, package = "MASS");
##' data=Cars93;
##' data=data[complete.cases(data),];
##' data=data[,sapply(data[1,],is.numeric)];
##' res=cont_selec(data,target="Horsepower",edge_cor=0.9,edge_explain=0.6,
##'                edge_param_number=3,family="gaussian",use_pforeach=FALSE);
##'
##' res_int=check_interaction(res,edge_cor=0.9,use_pforeach=FALSE)
##' plot_each_effect(res_int)
##' plot_combined_effect(res_int,sign_effect=-1)
##' plot_combined_effect(res_int,sign_effect=1)
##' 
##' 
#' @export
check_interaction <- function(res,edge_param_number=3,edge_cor=1.0,omit_contracted=TRUE,repeatn=50,use_pforeach=TRUE, perm=TRUE,check_spa_cor=FALSE){

    if(omit_contracted)res$result$sigmark[which(res$result$ngid>1)]="   ";
    result=res$result;
    
    data1=res$data1;
    data0=res$data0;
    family=res$family;
    group=res$group;
    ##if(flag_envs).ee.append("selec_stepwise",environment());
    list_sig=get_sig_list(res);
    vname=res$result$contname;
    ##data2=data1[,result$pv>-0.5];
    data2=data1;
    data02=data0;
    buf=data1[,1,drop=F];
    if(length(list_sig)>1){
        for(i in 1:(length(list_sig)-1)){
            for(j in (i+1):length(list_sig)){           
                buf=data1[,list_sig[i]]*data1[,list_sig[j]];
                data2=cbind(data2,buf);
                data02=cbind(data02,buf);
                myname=paste(vname[list_sig[i]],vname[list_sig[j]],sep="_XXX_");
                colnames(data2)[ncol(data2)]=myname;
                colnames(data02)[ncol(data02)]=myname;
                group$gid_data1=append(group$gid_data1,max(group$gid)+1);
                group$gid=append(group$gid,max(group$gid)+1);
                group$ngid=append(group$ngid,1);
                group$vname1_orig=append(group$vname1_orig,myname);
                group$vname0_orig=append(group$vname0_orig,myname);
                }
            }

    }
    print("number of explanatory variables:");
    print(ncol(data2)-1);
    print(colnames(data2));
    data2=scale_range(data2,minvalue=0.0,maxvalue=1.0);

 res2=selec_bestsub_stepwise(data02,data2,group, edge_param_number=edge_param_number,repeatn=repeatn,family=family, target=1, use_pforeach=use_pforeach,perm=perm,check_spa_cor=check_spa_cor,pos_x=res$pos_x,pos_y=res$pos_y);
res2$target=paste0(res$target,"_int");
  ## print(res2);

    check_cor(res2$data1[,res2$result$pv>-0.5],edge_cor=edge_cor);
              
    ##res=cont_selec(data2,pos_x=res$pos_x,pos_y=res$pos_y,target=1,edge_cor=edge_cor,edge_explain=edge_explain,edge_param_number=edge_param_number,repeatn=50,family=family,check_spa_cor=T,use_pforeach=TRUE,perm=TRUE);

    return(res2);
}



#' @export
sloppy_perm <- function(res0,nperm=5,check_order=NULL,seed=101,omit_pseudo=FALSE,outfile=NULL){
    if(flag_envs).ee.append("sloppy_perm",environment());
    
    if(length(outfile)==0){
        outfile="out_sloppy_perm";
    }
    
    pres=res0$result[,1:8];
    colnames(pres)[1]="sig";
    colnames(pres)[2]="contname";
    colnames(pres)[3]="n";
    colnames(pres)[4]="n2";
    colnames(pres)[5]="pv_slop";
    colnames(pres)[6]="pv_slop2";
    colnames(pres)[7]="pv_prr";
    colnames(pres)[8]="pv_paramet";
    
    pres$n=pres$n*0;
    pres$n2=pres$n2*0;    
    pres$pv_slop=pres$pv_slop*0;
    pres$pv_slop2=pres$pv_slop2*0;
    pres$pv_prr=res0$result$pv;
    pres$pv_paramet=res0$result$pv_paramet;

    
    list_sig=get_sig_list(res0);
    if(length(check_order)>0){
        if(sum(check_order==-1)>0){
            list_sig=list_sig[order(-1*seq(length(list_sig)))];
        }else{
                list_sig=list_sig[check_order];
            }

    }
    if(omit_pseudo){
        list_sig=list_sig[res0$result$ngid[list_sig]==1];
    }
    
    cat("#Check of explanatory variables with sloppy permutation#\n");
    cat("#Number of permutation: ",nperm,"#\n");
    cat("#Target variables:#\n");
    print(res0$result$contname[list_sig]);
    cat("\n");
    
    result1=cbind(res0$result,res0$result[,5:6]);    
    result1[,8:ncol(result1)]=res0$result[,6:ncol(res0$result)];
    result1[,6]=result1[,6]*0-1;
    result1[,7]=result1[,7]*0-1;
    colnames(result1)[6]="pv_slop";
    colnames(result1)[7]="pv_slop2";
    
    colnames(result1)[8:ncol(result1)]=colnames(res0$result)[6:ncol(res0$result)];
        
    
    
    for(k in 1:length(list_sig)){
        data1=res0$data1;
        permid=list_sig[k];
       
        expid=which(res0$result$pv>=0);
        expid=expid[-1*which(expid==permid)];
        
        if(length(expid)>0){
            expname=colnames(data1)[expid];
            depname=colnames(data1)[permid];
            form=reformulate(expname,depname);
            
            lm0=lm(form,data1);
            ##value=data1[,permid];
            var_perm_resi=lm0$residuals;
            coef=lm0$coefficients;
            var_perm_other=coef[1];
            for(i in 2:length(coef)){
                var_perm_other=var_perm_other+coef[i]*data1[,expid[i-1]];
            }
        }else{
            var_perm_resi=data1[,permid];
            var_perm_other=var_perm_resi*0.0;
            
        }

        
        data11=res0$data1[,res0$result$pv < -0.5];        
        data11=data11[,2:ncol(data11)];

        print("Checking regressor residuals:");
        cat(abs(data1[,permid]-(var_perm_other+var_perm_resi)),"\n");
        
        print("explanatory variable in the best model");
        print(colnames(res0$data1[,res0$result$pv > -0.5]));
        
        print("explanatory variable NOT in the best model");
        print(colnames(data11));
        
        set.seed(seed);
        for(i in 1:nperm){

            pres=sloppy_perm_sub(i,res0,permid,var_perm_resi,var_perm_other,data11,pres);
            cat("target:",res0$target,"\n");
            cat("#Trial: ",i," ",colnames(data1)[permid],"#\n");
            cat("# n:",pres$n[permid]," ",pres$n2[permid]," pv_slop:",pres$pv_slop[permid]," ",pres$pv_slop2[permid], "#\n");
  
            print(pres[append(1,list_sig),]);
          
            cat("\n perm exp var:",colnames(data1)[expid],"\n\n");


        }
        
        result1$pv_slop[permid]=pres$pv_slop[permid];
        result1$pv_slop2[permid]=pres$pv_slop2[permid];
    }

    result1$pv_slop[1]=nperm;
    result1$pv_slop2[1]=nperm;
    pres=pres[list_sig,];	
    pres=cbind(target=rep(res0$target,nrow(pres)),pres);

    outfilename=paste0(outfile,"_",res0$target,"_pv_slop.csv");
    write.table(result1,file=outfilename,sep=",",row.names=F);
    cat("output filename: ",outfilename,"\n");

    return(list(result=result1,pres=pres));
    
    
}

#' @export
sloppy_perm_sub <- function(i,res0,permid,var_perm_resi,var_perm_other,data11,pres){
            
    data1=res0$data1; 
    data0=res0$data0;
            

    var_perm_resi1=var_perm_resi;            
    var_perm_resi1=var_perm_resi1[sample(1:length(var_perm_resi1),replace=F)];

    ##buf0=var_perm_other+var_perm_resi;
    
    buf=var_perm_other+var_perm_resi1;
    buf=(buf-min(buf))/(max(buf)-min(buf));    
    
    cat("\n");
    print(buf);

    data1[,permid]=buf;
    
    permid_data0=which(colnames(data0)==res0$result$vname0_maxcor[permid]);
    data0[,permid_data0]=buf;
        
    pos_x=res0$pos_x;
    pos_y=res0$pos_y;
    group=res0$group;
    target=res0$target;
    
    res=selec_bestsub_stepwise(data0,data1,group,pos_x,pos_y,
                               edge_param_number=pa$edge_param_number,
                               repeatn=pa$repeatn,family=pa$family,
                               check_spa_cor=pa$check_spa_cor,target=target,
                               use_pforeach=pa$use_pforeach,perm=FALSE,nrep=pa$nrep);
            

    cond_step=as.integer((res$result$nmodel)*(res$result$stepwise)>0);
    cond_oneside=as.integer((res$result$nmodel)*(sign(res0$result$nmodel))==res$result$nmodel[1]);
    cond_twoside=as.integer(abs(res$result$nmodel)==res$result$nmodel[1]);
    
    condi=(cond_step*cond_oneside)[permid];
    condid=(cond_step*cond_twoside)[permid];
    
            
    if(condi+condid > 0){
        expname1=res$result$contname[-1*permid];
        expname1=expname1[res$result$pv[-1*permid]>=0];
        aic1drop=calc_aic(expname1,data=data1,family=res$bestmodel$family);
        
        expname00=res0$result$contname[-1*permid];
        expname00=expname00[res0$result$pv[-1*permid]>=0];        
        aic0drop=calc_aic(expname00,data=res0$data1,family=res0$bestmodel$family);

        aic1best=res$bestmodels$aic[1];
        aic0best=res0$bestmodels$aic[1];
        
        del_aic0=aic0drop-aic0best;
        del_aic1=aic1drop-aic1best;
        
        ##del_aic0=-1*res0$outer_bests$delta_aic[permid];
        ##del_aic1=-1*res$outer_bests$delta_aic[permid];
                
        passed=as.integer(del_aic1>del_aic0)+0.5*as.integer(del_aic1==del_aic0);
        
        if(condi>0)pres$n[permid]=pres$n[permid]+passed;
        if(condid>0)pres$n2[permid]=pres$n2[permid]+passed;
    }
            
    pres$n[1]=i;
    pres$n2[1]=i;
    
    pres$pv_slop[permid]=pres$n[permid]/as.double(i);
    pres$pv_slop2[permid]=pres$n2[permid]/as.double(i);
    
    pres$pv_prr[permid]=res0$result$pv[permid];
    
    cat("family is ",res$bestmodel$family,"\n");
    print(res$result[c(1,permid),2:4])
    
    if(condi+condid>0)cat("delta_aic0: ",del_aic0,"   delta_aic (resampled): ",del_aic1,"   difference: ",del_aic0 - del_aic1,"\n");
    if(condi+condid==0){
        cat("\n#condition not satisfied#\n");

    }
    
    return(pres);
}

#' @export
calc_aic <- function(expname,data,family){
    if(flag_envs).ee.append("calc_aic",environment());

    data_buf=cbind(data,id=seq(nrow(data)));   
    if(length(expname)>0)form = reformulate(expname, "y");
    if(length(expname)==0)form = as.formula("y ~ 1");
                  
    print(form);

    if(family=="glmm_poisson"){
        gg = glmmML::glmmML(formula=form,data=data_buf,family=poisson,cluster=id,start.sigma=0.001);        
    }
    
    if(family=="glmm_poisson_fixs"){
        gg = glmmML::glmmML(formula=form,data=data_buf,family=poisson,cluster=id,fix.sigma=TRUE,start.sigma=0.0);        
    }
    
    if((family=="gaussian")||(family=="poisson")){
        gg=glm(formula=form,data=data,family=family);
                
    }

    return(gg$aic);
}




#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/_/_/_/ Functions FOR PLOTTING _/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#


##' This list is a default parameter set for "plot_each_effect" and "plot_combined_effect" to plot result output of "cont_selec", "selec_bestsub_stepwise" and "check_interaction".
##' @title Default parameter set for "plot_each_effect" and "plot_combined_effect".
##' @author Hiroshi C. Ito
#' @export
param0=list(
    font_family="Times",
    ps=11,
    ncols=4,
    graph_edge_width=0.7,
    graph_bar_edge_width=0.2,
    graph_each_effect_xlab="Intensity",
    graph_each_effect_ylab="y",
    graph_each_effect_winh=3,
    graph_each_effect_winw=2.7,
    graph_each_effect_ymax_offset=1.2,
    graph_each_effect_bar_width=0.06,
    graph_each_effect_margin=data.frame(bottom=4.3,left=3.5,top=4.7,right=1.5),
    graph_combined_effect_xlab="Sample",
    graph_combined_effect_ylab="Combined effect",
    graph_combined_effect_vname_former=NULL,
    graph_combined_effect_vname_latter=NULL,
    graph_combined_effect_winh=3.5,
    graph_combined_effect_winw=4.5,
    graph_combined_effect_ymax_offset=1.5,
    graph_combined_effect_bar_width=0.4,    
    graph_combined_effect_margin=data.frame(bottom=3,left=4,top=4,right=2),
    fcol=NULL,
    name_swap=NULL
    );



##' This function plots separately the expected effect of each statistically-strong explanatory variable detected by "cont_selec" or "selec_bestsub_stepwise". See "cont_selec" in this manual or help(cont_selec).
##' @title plots effects of statistically-strong explanatory variables separately
##' @param res list : result output of "cont_selec"  or "selec_bestsub_stepwise"
##' @param param list : graphics parameter. Default values are given by "param0"
##' @param outfile character : filename of image output.
##' @param outpng TRUE/FALSE : if TRUE, plottings are saved in png format.
##' @param outeps TRUE/FALSE : if TRUE, plottings are saved in eps format.
##' @param outpng_show  TRUE/FALSE : if TRUE, saved png files are shown.
##' @param hide_title  TRUE/FALSE : if TRUE, titles of plottings are not shown.
##' @param sign_flip  TRUE/FALSE : if TRUE, signs of mean and max impacts of explanatory variables are flipped.
##' @param pch  integer or "text" : pch for schatter plot. If "text", then sample ids are plotted.
##' @author Hiroshi C. Ito
#'@export
plot_each_effect <- function(res,param=param0,outfile="testout",outpng=F,outeps=F,outpng_show=T,hide_title=F,sign_flip=FALSE,pch=24){

    if(flag_envs).ee.append("plot_each_effect",environment());

    if(length(get_sig_list(res))==0){
        print("no statistically-strong variable!!");
        return();
    }
    family=res$family;
    result=res$result;
    sigmark=result$sigmark;
    lists=which(sigmark!="   ");
    listn=which(as.logical((sigmark=="   ")*(result$pv>0)));    

    data1=res$data1;
    coef=result$coef;
    coef_se=result$coef_se;
    names(coef)=result$contname;
    vname=result$contname;
    vname1=adjust_vname(result,param);
    fcol=make_fcol(param$fcol,vname);

    c=calc_plot_each(coef,coef_se,data1,listn,lists,param,family,sign_flip);
   
    offset=c$offset;
    xxx=c$xxx;
    yyy=c$yyy;
    yest0=c$yest0;
    yestm=c$yestm;
    yest1=c$yest1;
    maximp=c$maximp;
    meanimp=c$meanimp;
        
    xcurv=c$xcurv;
    ycurv=c$ycurv;
    ylim=c$ylim;
    
    ymax=max(yyy);
    ymin=min(yyy);

    dev_offset= cur.dev();
    
    for(i in 1:length(lists)){
        col.main=fcol[lists[i]];
        ##text_impact=sprintf("max: %2.2f (%s)  mean: %2.2f",maximp[i],sigmark[lists[i]],meanimp[i]);
        text_impact=sprintf("max : %2.2f (%2.2f ~ %2.2f)\nmean: %2.2f (%2.2f ~ %2.2f)",maximp[i],c$maximp_lower[i],c$maximp_upper[i],meanimp[i],c$meanimp_lower[i],c$meanimp_upper[i]);
        ##text_impact1=sprintf("max: %2.2f (%2.2f ~ %2.2f)",maximp[i],c$maximp_lower[i],c$maximp_upper[i]);
        ##text_impact2=sprintf("mean: %2.2f (%2.2f ~ %2.2f)",meanimp[i],c$meanimp_lower[i],c$meanimp_upper[i]);
        
   
        text_title=sprintf("%s\n%s",res$target,vname1[lists[i]]);
        if(hide_title)text_title=sprintf("%s",vname1[lists[i]]);
            
        flag_expr=0;
        if(length(grep("_XXX_",vname1[lists[i]]))>0){
            flag_expr=1;
            titles=strsplit(vname1[lists[i]],split="_XXX_")[[1]];
            ##text_title= paste("expression(atop(",res$target,",",titles[1],"%*%",titles[2],"))");
            text_title= sprintf("expression(atop(bold(%s), bold(\"%s\" %%*%% \"%s\")))",res$target,titles[1],titles[2]);
            if(hide_title)text_title=sprintf("expression(bold(\"%s\" %%*%% \"%s\"))",titles[1],titles[2]);
            
            text_title=eval(parse(text=text_title));
        }

        set_plot_each(param,ylim,text_title,col.main,flag_expr);
                
        x=c(0,mean(xxx[,i]),1);
        ytop=c(yest0[i],yestm[i],yest1[i]);
        barw=param$graph_each_effect_bar_width;
        ##        colmin=rgb(0.89,0.9,0.9);
        ##        colmean=rgb(0.49,0.5,0.5);
        ##        colmax=rgb(0.0,0.01,0.01);
        colmin=rgb(1,1,1);
        colmean=rgb(0.89,0.9,0.9);
        colmax=rgb(0.49,0.5,0.5);
        
        rect(xleft=x-barw, ybottom=0, xright=x+barw, ytop=ytop, col=c(colmin,colmean,colmax));
        
     
        if(pch!="text"){
            points(xxx[,i],yyy[,i],pch=pch,cex=0.7,bg="white",lwd=0.5);        
        }

        lines(xcurv[,i],ycurv[,i],col=col.main,lwd=2);

       if(pch=="text"){
           text(x=xxx[,i],y=yyy[,i],labels=as.character(1:length(xxx[,i])),cex=0.7,col="black",font=1);
##           text(x=xxx[,i],y=yyy[,i],labels=as.character(1:length(xxx[,i])),cex=0.7,col=col.main,font=1);
       }
       ##mtext(text_impact, side = 3, line = 0.1, adj=1,at = NA);
        ##if(hide_title)
            mtext(text_impact, side = 3, line = 0.1, adj=0,at = NA);
       ## if(!hide_title){
       ##     mtext(text_impact1, side = 3, line = 0.1, adj=0,at = NA);
       ##     mtext(text_impact2, side = 4, line = 0.1, adj=0,at = NA);
       ##     }
        tickValues = c(0,0.5,1)
        tickStrings = sprintf("%.1f",tickValues) 
        axis(side=1,at=tickValues,labels=tickStrings,cex=1.0)
        box(lty=1);
    }

    vname=res$result$contname;
    if(outpng||outeps){
        for(i in 1:length(lists)){
            ofilename=paste(outfile,"_",res$target,"_",vname[lists[i]],sep="");
            pngout(dev_offset+i,ofilename,outeps=outeps,show=outpng_show,font=param$font_family,outpng=outpng);
        }
    }
}

#'@export
get_sig_list <- function(res){
    return(which(res$result$sigmark!="   "));
}


calc_plot_each <- function(coef,coef_se,data1,listn,lists,param,family,sign_flip){
    
    resi=get_resi0(coef,data1,family=family)                 

    offset=coef[1];
    if(length(listn)>0){
        for(i in 1:length(listn)){
            offset=offset+coef[listn[i]]*mean(data1[,listn[i]]);
        }
    }
  
    
    nn=length(lists);
    
    yest0=yestm=yest1=maximp=meanimp=array(0.0,nn);
    maximp_upper=maximp_lower=maximp;
    meanimp_upper=meanimp_lower=meanimp;
    
    xxx=matrix(0.0,nrow=nrow(data1),ncol=nn);
    
    yyy=xxx;
    xcurv=matrix(0.0,nrow=64,ncol=nn);
    ycurv=xcurv;
    for(i in 1:length(lists)){
        mycoef=coef[lists[i]];
        myse=coef_se[lists[i]];
        xcurv[,i]=seq(0,1,length=64);
        ycurv[,i]=(offset+mycoef*xcurv[,i]);
        
        xx=data1[,lists[i]];
        yest=(offset+mycoef*xx);
        yest0[i]=(offset+mycoef*0);
        yestm[i]=(offset+mycoef*mean(xx));
        yest1[i]=(offset+mycoef*1);

        if(sign_flip){
            mycoef= -1*mycoef;
            }

        maximp[i]=mycoef;
        meanimp[i]=mycoef*mean(xx);

        maximp_upper[i]=mycoef+1.96*myse;
        maximp_lower[i]=mycoef-1.96*myse;
        meanimp_upper[i]=(mycoef+1.96*myse)*mean(xx);
        meanimp_lower[i]=(mycoef-1.96*myse)*mean(xx);
        
        meanimp[i]=mycoef*mean(xx);
     
        
        if((family=="glmm_poisson")||(family=="poisson")){
            ycurv[,i]=exp(ycurv[,i]);
            yest=exp(yest);
            yest0[i]=exp(yest0[i]);
            yestm[i]=exp(yestm[i]);
            yest1[i]=exp(yest1[i]);
            
            maximp[i]=exp(maximp[i]);
            meanimp[i]=exp(meanimp[i]);
            maximp_upper[i]=exp(maximp_upper[i]);
            maximp_lower[i]=exp(maximp_lower[i]);
            meanimp_upper[i]=exp(meanimp_upper[i]);
            meanimp_lower[i]=exp(meanimp_lower[i]);
            
            
        }
        yyy[,i]=yest+resi;
        xxx[,i]=xx;
    }


    ymin=min(yyy);
    ymax=max(yyy);
    if((family=="glmm_poisson")||(family=="poisson"))ylim=c(0,ymax*param$graph_each_effect_ymax_offset);
    if(family=="gaussian")ylim=c(min(ymin,ymax),max(ymin,ymax))*param$graph_each_effect_ymax_offset;
    
    return(list(offset=offset,xxx=xxx,yyy=yyy,yest0=yest0,yestm=yestm,yest1=yest1,maximp=maximp,meanimp=meanimp,xcurv=xcurv,ycurv=ycurv,ylim=ylim,resi=resi,maximp_upper=maximp_upper,maximp_lower=maximp_lower,meanimp_upper=meanimp_upper,meanimp_lower=meanimp_lower));
}



set_plot_each <- function(param,ylim,text_title,col.main,flag_expr){
    if(flag_X11){
        X11(width=param$graph_each_effect_winw,height=param$graph_each_effect_winh,bg="white");
    }else{
        dev.new(width=param$graph_each_effect_winw,height=param$graph_each_effect_winh,bg="white");
    }
        
        par(mar=as.numeric(param$graph_each_effect_margin),ps=param$ps)
        par(family=param$font_family);
        par(bty="o", yaxs="i",mgp=c(2,0.7,0)) ;
        
        xlab_each=param$graph_each_effect_xlab;
    ylab_each=param$graph_each_effect_ylab;

    
    plot(c(1,1),c(1,1),type="n",xlab=xlab_each,cex.lab=1.2,ylab=ylab_each,xaxt="n",xlim=c(-0.1,1.1),ylim=ylim,cex.main=1.2,cex.axis = 1.0);
    
    ##mtext(text_title, side = 3, line = 1.0, adj=0.5,at = NA,cex=1.2,col=col.main,font=2);
    ##mtext(text_title, side = 3, line = 1.1-0.2*flag_expr, adj=0.5,at = NA,cex=(1.2),col=col.main,font=2);
    
    mtext(text_title, side = 3, line = 2.0-0.2*flag_expr, adj=0.5,at = NA,cex=(1.2),col=col.main,font=2);
    rect(xleft=-1.0,ybottom=-1,xright=2.0,ytop=ylim[2]*10,col="white");
}

adjust_vname <- function(result,param){
        vname=result$contname;
        mask=result$ngid>2;
        vname[mask]=paste0(vname[mask],"(",result$vname0_maxcor[mask],")");
        mask=result$ngid==2;
        vname[mask]=result$vname0_maxcor[mask];
        vname1=swap_vname(vname,param$name_swap);
        return(vname1)
}

swap_vname <- function(vname,name_swap){
    vname1=vname;
    if(!is.null(name_swap)){
        for(i in 1:length(name_swap)){
            vname1[vname==names(name_swap)[i]]=name_swap[[i]];
        }
    }
    return(vname1);
}

make_fcol <- function(fcol0=NULL,vname,flip=FALSE){
    if(is.null(fcol0)){
        fcol=rep("",length(vname));
        h=0.0;
        s=1.0;
        v=0.9;
        ds=0.7*5.0/length(vname);
        for(i in 1:length(fcol)){
            fcol[i]=hsv(h,s,v);
            h=h+0.15;
            if(h>1.0)h=h-1.0;
            if(i%%6==0){
                s=max(s-ds,0.2);
                v=max(v-ds*0.8,0.2);
            }
        }
        if(flip)fcol=fcol[length(fcol):1];
        
    }else{        
            fcol=rep("white",length(vname));
            for(i in 1:length(fcol0)){
                fcol[vname==names(fcol0)[i]]=fcol0[[i]];
            }
    
        }
        return(fcol);
}


cur.dev <- function(){
    v=as.numeric(dev.list()[length(dev.list())]);
    if(length(v)==0)v=1;
return(v);
}


##' This function plots the expected combined effect of statistically-strong explanatory variables detected by "cont_selec" or "selec_bestsub_stepwise". See "cont_selec" in this manual or help(cont_selec).
##' @title plots combined effect of statistically-strong explanatory variables
##' @param res list : result output of "cont_selec" or "selec_bestsub_stepwise"
##' @param param list : graphics parameter. Default values are given by "param0"
##' @param outfile character : filename of image output.
##' @param outpng TRUE/FALSE : if TRUE, plottings are saved in png format.
##' @param outeps  TRUE/FALSE : if TRUE, plottings are saved in eps format.
##' @param outpng_show  TRUE/FALSE : if TRUE, saved png files are shown.
##' @param hide_title  TRUE/FALSE : if TRUE, titles of plottings are hidden.
##' @param sign_effect  1/-1 : if 1, only explanatory variables with positive regression coeffecients are shown, and vice versa.
##' @param do_order  TRUE/FALSE : if TRUE, ponds are ordered along horisontal axis in terms of maximum combined impacts.
##' @author Hiroshi C. Ito
#' @export                                
plot_combined_effect <- function(res,param=param0,outfile="testout",outpng=F,outeps=F,outpng_show=T,hide_title=F,sign_effect=1,do_order=TRUE){

    if(flag_envs).ee.append("plot_combined_effect",environment());
    
    if(length(get_sig_list(res))==0){
        print("no statistically-strong variable!!");
        return();
    }

    family=res$family;
    result=res$result;
    sigmark=result$sigmark;
    
    data1=res$data1;
    coef=result$coef;
    coef_se=result$coef_se;
    names(coef)=result$contname;
    vname=result$contname;
    vname1=adjust_vname(result,param);
    fcol=make_fcol(param$fcol,vname);
    lists=which((sigmark!="   ")*coef*sign_effect>0);
          
    resi=get_resi0(coef,data1,family=res$family)                 
    r2=1.0-var(resi)/var(data1[,1]);

    listss=adjust_lists_order(lists,fcol,vname,param);
    ffcols=fcol[listss,drop=F];
  
    mmm=data1[,listss,drop=F];
    
    mmm_se=data1[,listss,drop=F];
    data1_mean=colMeans(mmm);
    mmm_mean=mmm[1,];
    mmm_mean_se=mmm[1,];
    for(j in 1:ncol(mmm)){
        mmm[,j]=sign_effect*coef[listss[j]]*data1[,listss[j]];
        mmm_se[,j]=(coef_se[listss[j]]*data1[,listss[j]])^2;
        mmm_mean[j]=sign_effect*coef[listss[j]]*data1_mean[j];
        mmm_mean_se[j]=(coef_se[listss[j]]*data1_mean[j])^2;
    }

  
    impa=(apply(mmm,1,cumsum));
    ymax_id=which.max(rowSums(mmm));
    ymax=(rowSums(mmm))[ymax_id];
    ymax_se=sqrt(rowSums(mmm_se))[ymax_id];
    ymax_upper=ymax+1.96*ymax_se;
    ymax_lower=ymax-1.96*ymax_se;
    
    
    
    ymean=sum(mmm_mean);
    ymean_se=sqrt(sum(mmm_mean_se));
    ymean_upper=ymean+1.96*ymean_se;
    ymean_lower=ymean-1.96*ymean_se;
    
    impa1=rbind(rep(0,nrow(data1)),impa);
    
    if((family=="glmm_poisson")||(family=="poisson")){
        impa=exp(impa);
        ymax=exp(ymax);
        ymax_upper=exp(ymax_upper);
        ymax_lower=exp(ymax_lower);
        ymean=exp(ymean);
        ymean_upper=exp(ymean_upper);
        ymean_lower=exp(ymean_lower);

        impa1=exp(impa1);
    }
     
    if(length(listss)>1)yymax=impa[length(listss),];
    if(length(listss)==1)yymax=impa;
    
    nsamp=nrow(data1);
    nvar=(sum(res$result$nmodel!=0)-1);
    set_plot_combined(param,res$target,r2,nvar,ymax,ymax_upper,ymax_lower,ymean,ymean_upper,ymean_lower,nsamp,family,sign_effect,hide_title);

    barw=param$graph_combined_effect_bar_width;;    

    oo=order(-1*yymax);
    if(do_order==FALSE)oo=1:length(oo);
    x=1:nrow(data1);
   
    for(i in 1:length(listss)){
        yb=impa1[i,];
        yt=impa1[(i+1),];
        rect(xleft=x-barw, ybottom=yb[oo], xright=x+barw, ytop=yt[oo], col=ffcols[i],lwd=param$graph_bar_edge_width);
    }

    text_legend=vname1[listss[length(listss):1]];
    for(i in 1:length(text_legend)){
        if(length(grep("_XXX_",text_legend[i]))>0){
            titles=strsplit(text_legend[i],split="_XXX_")[[1]];
            text_buf=text_title=sprintf("expression(bold(\"%s\" %%*%% \"%s\"))",titles[1],titles[2]);
            text_legend[i]=eval(parse(text=text_buf));
        }
    }
    
    legend("topright", legend = text_legend, col = "black",fill=ffcols[length(ffcols):1],cex=0.9);
    ikeid=as.character(1:nrow(data1));
    text(x,yymax[oo],label=ikeid[oo],pos=3,cex=0.8);    
    
    box(lty=1);
    
    if(outpng||outeps){
        ofilename=paste(outfile,"_",res$target,"_combined",sep="");
        pngout(cur.dev(),ofilename,outeps=outeps,show=outpng_show,outpng=outpng,font=param$font_family);        
    }
}

adjust_lists_order <- function(lists,fcol,vname,param){
    vname_former=param$graph_combined_effect_vname_former;
    vname_latter=param$graph_combined_effect_vname_latter;
    
    vnameb=vname[lists];
    mask_vnameb=1:length(vnameb);
       
    if(!is.null(vname_former)){
        for(i in 1:length(vname_former))mask_vnameb[vnameb==vname_former[i]]= -1*i;
    }
    if(!is.null(vname_latter)){
        for(i in 1:length(vname_latter))mask_vnameb[vnameb==vname_latter[i]]=(200+i); 
    }
    listss=lists[order(mask_vnameb)];
    return(listss);
}


set_plot_combined <- function(param,target,r2,nvar,ymax,ymax_upper,ymax_lower,ymean,ymean_upper,ymean_lower,nsamp,family,sign_effect,hide_title){

    text_title=sprintf("%s (%2.1f%% explaind by %dvars, %s)\n\n",target,100*r2,nvar,family);
    if(hide_title)text_title=sprintf(" \n\n");
    text_impact=sprintf("max: %2.2f (%2.2f ~ %2.2f)  mean: %2.2f (%2.2f ~ %2.2f)",ymax,ymax_lower,ymax_upper,ymean,ymean_lower,ymean_upper);
    
    ylab_combined=param$graph_combined_effect_ylab;
    xlab_combined=param$graph_combined_effect_xlab;
    if(flag_X11){
        X11(width=param$graph_combined_effect_winw,height=param$graph_combined_effect_winh);
    }else{
        dev.new(width=param$graph_combined_effect_winw,height=param$graph_combined_effect_winh);
    }
    par(mar=as.numeric(param$graph_combined_effect_margin),ps=param$ps)
    par(lwd=param$graph_edge_width);
    par(family=param$font_family);
    
    par(yaxs="i",mgp=c(1.8,0.7,0));
    
    if(sign_effect<0)ylab=paste0(ylab_combined," (negative)");
    if(sign_effect>0)ylab=paste0(ylab_combined," (positive)");
    
    if((family=="glmm_poisson")||(family=="poisson")){
        ylim1=ymax^(param$graph_combined_effect_ymax_offset);
        plot(1:nsamp,type="n",log="y",ylim=c(1,ylim1),ylab=ylab_combined,xlab="",xaxt="n",cex.axis = 0.9,main=text_title,cex.main=0.9,cex.lab=1.2);
    }
    if(family=="gaussian"){
        ylim1=ymax*(param$graph_combined_effect_ymax_offset);
        plot(1:nsamp,type="n",ylim=c(0,ylim1),ylab=ylab,xlab="",xaxt="n",cex.axis = 0.9,main=text_title,cex.main=0.9,cex.lab=1.2);
    }
    
    mtext(text_impact, side = 3, line = 0.1, adj=0,at = NA);
    mtext(xlab_combined, side = 1, line = 0.5, adj=0.5,at = NA,cex=1.2);
}



##' This function saves plotting in png and/or eps format by using "dev.copy2eps" and command "convert" in ImageMagick.
##' @title saves plotting in png format
##' @param dev_id integer : device id (default value is current device id)
##' @param plotfile character : filename
##' @param density integer : "density" option for "convert"
##' @param geometry integer : "geometry" option for "convert"
##' @param outeps TRUE/FALSE : if TRUE plotting is saved in eps format
##' @param prefix character : location of saved file
##' @param show TRUE/FALSE : if TRUE saved png file is shown by "display" command.
##' @param font character : font family
##' @param outpng  TRUE/FALSE : if TRUE plotting is saved in eps format
##' @author Hiroshi C. Ito
#' @export
pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=100*100.0/density,outeps=FALSE,prefix="./",show=TRUE,font=param0$font_family,outpng=TRUE){
    
    dens=as.character(density);
    geom=as.character(as.integer(geometry));
    dev.set(dev_id);

    dev.copy2eps(file=".R_pngout_temp.eps",family=font);
   
     if(outeps){
        system(paste("cp .R_pngout_temp.eps ",prefix,plotfile,".eps",sep=""));
        cat("eps output:",paste(prefix,plotfile,".eps",sep=""),"\n")
     }
    
    if(outpng){
        system(paste("convert -density ",dens,"x",dens," -geometry ", geom,"%","  -background white -alpha remove .R_pngout_temp.eps .R_pngout_temp.png",sep=""));
        system(paste("mv .R_pngout_temp.png ",prefix,plotfile,".png",sep=""));
        cat("png output:",paste(prefix,plotfile,".png",sep=""),"\n");

        if(show)system(paste("display ",prefix,plotfile,".png&",sep=""));
    }
    system("rm .R_pngout_temp.eps ");
    
}


