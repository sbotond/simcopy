##	
## Copyright 2009 Botond Sipos	
## See the package description for licensing information.	

##
## Constructor: SCRootSummary
##
##########################################################################/** 
#
# @RdocClass SCRootSummary
# 
# @title "The SCRootSummary class"
# 
# \description{ 
#	SCRootSummary objects are blessed lists containing summary entries created by
#	\code{summary.*} methods.
#	
#	@classhierarchy
# }
#	
# @synopsis
#	
# \arguments{
# 	\item{summary}{A list.}
# 	\item{...}{Not used.}
#	}
# 
# \section{Fields and Methods}{ 
# 	@allmethods
# }
# 
# @author
#
# \seealso{ 
# 	@seeclass 
# }
# 
#*/###########################################################################
setConstructorS3(
  "SCRootSummary",
  function(summary=list(),...){
			
			# Stepping out of the R.oo framework to provide 
			# the expected behaviour.
			class(summary)<-c("SCRootSummary");
			summary;
  },
  ###
  enforceRCC=FALSE
);

##
## Method: print.SCRootSummary
##
###########################################################################/**
#
# @RdocMethod print
# 
# @title "Print out a SCRootSummary object" 
# 
# \description{ 
#	@get "title".
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{x}{A SCRootSummary object.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	The summary object (invisible).
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "print",
  class="SCRootSummary",
  appendVarArgs=FALSE,
  function(
    x,
    ...
  ){
	this<-x;
	cat("\n");
	for (i in names(this)){
        cat(paste(i,": ",this[[i]],"\n",sep=""));
   	}
		cat("\n");
		invisible(this);

  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##########################################################################/** 
#
# @RdocClass SCRoot
# 
# @title "The root class for all simcopy objects"
# 
# \description{ 
#		The root class for all simcopy objects containig some utility methods. 
#		@classhierarchy
# }
#	
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#	
# \section{Fields and Methods}{ 
# 	@allmethods
# }
# 
# \examples{ 
#		obj<-SCRoot();
#		ll(obj);
# }
# 
# @author
#
#
# \seealso{ 
# 	Object
# }
# 
#*/###########################################################################
setConstructorS3(
  "SCRoot",
  function(...){
  extend(Object(), "SCRoot",
		.comments=character(0),
		.summary=list()
  );
  },
  ###
  enforceRCC=TRUE
);

##
## Method: getMethodsList
###########################################################################/**
#
# @RdocMethod getMethodsList
# 
# @title "Get the list of applicable methods for an object" 
# 
# \description{ 
#	@get "title".
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{A SCRoot object.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	The list of applicable methods.
# } 
# 
# \examples{
#	# create an object
#	o<-SCRoot()
#	# get the applicable methods
#	getMethodsList(o)
#	# get methods via virtual field
#	o$methodsList
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
##
setMethodS3(
  "getMethodsList",
  class="SCRoot",
  function(
    this,
    ...
  ){

			class<-class(this)[[1]];
			mlist<-getMethods.Class(this);

			# If the class has no methods, do not 
			# consider the methods from the parent class.
			if(names(mlist)[[1]] == class){	
      			as.character(names(mlist[[1]]));
			}
			else {
				return(character(0));
			}

  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: ll
##
###########################################################################/**
#
# @RdocMethod ll
# 
# @title "Display detailed information about the virtual fields and methods defined for a given object" 
# 
# \description{ 
#	@get "title".
#	The method prints the class of the object, all the parent classes,
#	and the virtual fields and methods defined in the immediate class.
#
#	This method provides a "quick and minimal" documentation for SimCopy classes.
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{A SCRoot object.} 
# 	\item{quiet}{Do not print out methods list.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	Text.
# } 
# 
# \examples{
#	# create a Site object
#	s<-Site()
#	ll(s)
#	# get information about the Process class
#	ll(Process())
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "ll",
  class="SCRoot",
  function(
    this,
    quiet=FALSE,
    ...
  ){
		
		class<-class(this);
		parents<-class[-1];
		class<-class[[1]]
		methods<-getMethodsList(this);
		fields<-getFields(this);
		text<-character(0);	

		pretty.print<-function(vec,text){

				tmp<-"";
				if(length(vec) > 0 ){
					tmp<-paste(tmp,"  ",vec,sep="",collapse="\n");
				  paste(text,tmp,"\n",sep="");
				} else {
					return(text);
				}
		}

	
		text<-paste(text,"\nClass: ",class,"\n",sep="");
		text<-paste(text,"Inherits from: ",paste(parents,collapse=" "),"\n",sep="");
		text<-paste(text,"Fields (",length(fields),"):\n",sep="");
		text<-pretty.print(fields,text);	

		# Discriminate between the methods implementing 
		# virtual fileds and the rest:
	
		vfields<-character(0);
		methods.not.virtual<-character(0);

		num.args<-function(fun){
			length(formals(fun))
		}

		method.to.field<-function(method){

			 method<-sub('^(get|set)(.*)','\\2',method);
			 tmp<-as.array(strsplit(method,"",fixed=TRUE))[[1]];
       tmp[1]<-tolower(tmp[1]);
       paste(tmp,collapse="");			

		}

		classify.method<-function(method,limit) {

				if( num.args( paste(method,".",class(this)[[1]],sep="") ) == limit){
                vfields<<-c(vfields,method.to.field(method));
            } else {
              methods.not.virtual<<-c(methods.not.virtual,method);
            }

		}

		for(method in methods){
			
				# Get methods for virtual fields have 2 aguments: "this" and "...".
				if(length(grep("^get",method,perl=TRUE)) == 1) {
					classify.method (method,limit=2)
				}
				# Set methods for virtual fields have 3 aguments: "this", "..." and "value".
				else if (length(grep("^set",method,perl=TRUE)) == 1) {
					classify.method (method,limit=3)
				} else {
					methods.not.virtual<-c(methods.not.virtual,method);
				}
		
		}
		vfields<-sort(unique(vfields));	

		lapply(methods.not.virtual,
			function(name) {
				tmp<-method.to.field(name);
				if (length(intersect(tmp,vfields)) > 0 ) {
					print(intersect(tmp,vfields));
					throw("Method classification inconsistency! Blaming ",paste(intersect(tmp,vfields),collapse=" "),". \n");
				}
			}
		);
		
		text<-paste(text,"Virtual fields (",length(vfields),"):\n",sep="");
		text<-pretty.print(vfields,text);
		text<-paste(text,"Methods implemented in ",class," (",length(methods.not.virtual),"):\n",sep="");
		text<-pretty.print(sort(methods.not.virtual),text);
		text<-paste(text,"\n",sep="");
		
		if(!quiet){ cat(text) }	

		invisible(text);

  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: getComments
##
###########################################################################/**
#
# @RdocMethod getComments
# 
# @title "Get the comments associated with an object" 
# 
# \description{ 
#	@get "title".
#
#	The comment field can contain any type of object.
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{A SCRoot object.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	The value of the comment field.
# } 
# 
# \examples{
#	# create an object
#	o<-SCRoot()
#	# add some comments
#	setComments(o,"Random comment")
#	# get the comment 
#	getComments(o)
#	# get/set the comment via virtual fiels
#	o$comments<-"Second random comment"
#	o$comments
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "getComments",
  class="SCRoot",
  function(
    this,
    ...
  ){
			this$.comments;
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: setComments
##
###########################################################################/**
#
# @RdocMethod setComments
# 
# @title "Set the comments associated with an object" 
# 
# \description{ 
#	@get "title".
#
#	The comment field can contain any type of object.
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{A SCRoot object.} 
#	\item{new_value}{An object.}
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	The new value of the comment field (invisible).
# } 
# 
# \examples{
#	# create an object
#	o<-SCRoot()
#	# add some comments
#	setComments(o,"Random comment")
#	# get the comment 
#	getComments(o)
#	# get/set the comment via virtual fiels
#	o$comments<-"Second random comment"
#	o$comments
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "setComments",
  class="SCRoot",
  function(
    this,
    new_value,
    ...
  ){
			this$.comments<-new_value;
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: summary.SCRoot
##
###########################################################################/**
#
# @RdocMethod summary
#
# @title "Summarize the properties of an object"
#
# \description{
#       @get "title".
# }
#
# @synopsis
#
# \arguments{
#       \item{object}{An object}
#       \item{...}{Not used.}
# }
#
# \value{
#  	Returns a SCRootSummary object.
# }
#
# \examples{
#
#       # create an object
#       a<-SCRoot()
#       # get a summary
#       summary(a)
# }
#
# @author
#
# \seealso{
#       @seeclass
# }
#
#*/###########################################################################
setMethodS3(
  "summary",
  class="SCRoot",
  function(
    object,
    ...
  ){
		this<-object;	
		# Adding the Comments field:
		if(length(this$.comments) > 0 ) {
		this$.summary$Comments<-paste(this$.comments, collapse=", ");
		}
		
		obj<-SCRootSummary(summary=this$.summary);
		this$.summary<-list();
		# Return a summary object:
		return(obj);

	},
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: is.na.SCRoot
##
###########################################################################/**
#
# @RdocMethod is.na
# 
# @title "Check if a SCRoot object is NA" 
# 
# \description{ 
#	@get "title".
#	SCRoot objects accanot be NA, so this method always returns FALSE.
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{x}{A SCRoot object.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	FALSE
# } 
# 
# \examples{
#	is.na(SCRoot());
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "is.na",
  class="SCRoot",
  function(
    x,
    ...
  ){
		
		# We don't want our objects to be NA-s!	
		return(FALSE);
		
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: plot.SCRoot
##

###########################################################################/**
#
# @RdocMethod plot
# 
# @title "Dummy plot method for the SCRoot class" 
# 
# \description{ 
#	@get "title".
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# FALSE
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "plot",
  class="SCRoot",
  function(
    ...
  ){
	
		cat("No plot method defined for this object!\n");	
		return(invisible(FALSE));
		
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

##
## Method: is.SCRoot.default
##
###########################################################################/**
#
# @RdocDefault is.SCRoot
# 
# @title "Check if an object inherits from SCRoot" 
# 
# \description{ 
#	@get "title".
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{An object.} 
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 	TRUE or FALSE.
# } 
# 
# \examples{
#	# create some objects
#	o<-SCRoot()
#	a<-Alphabet()
#	x<-Object()
#	# check if they inherit form SCRoot
#	is.SCRoot(o)
#	is.SCRoot(a)
#	is.SCRoot(x)
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "is.SCRoot",
  class="default",
  function(
    this,
    ...
  ){

		if(!is.object(this)) {return(FALSE)}
		inherits(this,"SCRoot");
		
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);
