#define maxdomain 100 //the maximum number of domains followed in class
							//substructure

class substructure {// this is a stack that keeps track of domains yet
							//	to be assigned coordinates in void coordinates
private:
int nuc5[maxdomain],num;
//nuc5 is the number of the 5' base leading into a domain
//num is the number of domains being stored

public:
bool getsub(int *getnuc5);
void putsub(int putnuc5);
substructure();
};

substructure::substructure(){
	num=0;
}

//	Put a 5' nuc and its x and y coordinates onto the stack
void substructure::putsub(int putnuc5){


	//cout << "putting "<<putnuc5<<"\n";

   nuc5[num]=putnuc5;
	num++;

   //cout << "num = "<<num<<"\n";
}

//	Return the next 5' nuc and its x and y coordinates
bool substructure::getsub(int *getnuc5){
  	if(num==0) return false;
   else{
   	num--;
      *getnuc5=nuc5[num];
      //cout << "sending "<<*getnuc5<<"\n";
      //cout << "num = "<<num<<"\n";
      return true;
   }
}
