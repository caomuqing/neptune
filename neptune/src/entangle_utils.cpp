/* ----------------------------------------------------------------------------
 * Copyright 2021, Cao Muqing, Internet of Things Laboratory
 * Nanyang Technological University
 * All Rights Reserved
 * Authors: Cao Muqing
 * -------------------------------------------------------------------------- */

#include "entangle_utils.hpp"
// #include <iostream>
#include "termcolor.hpp"
#include "gjk.hpp"

//this function indicates if ac is clockwise or anticlockwise to ab
// >0: anti-clockwise
// <0: clockwise
double eu::vectorWedge2(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
  return (b(0)-a(0))*(c(1)-a(1))-(c(0)-a(0))*(b(1)-a(1));
}

double eu::vectorWedge2(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, 
                        Eigen::Vector2d& ab, Eigen::Vector2d& ac)
{
  ab = b - a;
  ac = c - a;
  return ab(0)*ac(1)-ac(0)*ab(1);
}

// check if two lines (a1, a2) and (b1, b2) intersect within their segment
bool eu::lineSegIntersect(Eigen::Vector2d& a1, Eigen::Vector2d& a2, 
                          Eigen::Vector2d& b1, Eigen::Vector2d& b2)
{
  double c1 = vectorWedge2(b1, a1, a2);
  double c2 = vectorWedge2(b2, a1, a2);

  if (c1*c2<0)
  {
    double d1 = vectorWedge2(a1, b1, b2);
    double d2 = vectorWedge2(a2, b1, b2);    

    if (d1 * d2 <0) return true;
    else return false;
  }
  else return false;
}

// bool eu::checkEntangleUpdate(int& state, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik,
//                                Eigen::Vector2d& pikplus1, Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength,
//                                double maxheight)
// {
//   // std::cout<<termcolor::bold<<termcolor::green<<"checking entangling condition!"<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<"pk is"<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<pk<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<"pkplus1 is"<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<pkplus1<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<"pik is"<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<pik<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
//   // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";

//   double c1 = vectorWedge2(pk, pik, pbi);
//   double c2 = vectorWedge2(pkplus1, pikplus1, pbi);
//   double d1 = vectorWedge2(pk, pik, pb);
//   double d2 = vectorWedge2(pkplus1, pikplus1, pb);
//   if ((c1*c2)<0)
//   {

//     if (d1>0 && d2>0) //crossing in between agent i and its base
//     {
//       // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and its base!"<<"\n";
//       if (state==0) state = 1;
//       else if (state==1) state=0;
//       else return true; //entangle happens

//     }
//     else if (d1<0 && d2<0) //crossing in between agent i and the workspace boundary
//     {
//       // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and workspace boundary!"<<"\n";
      
//       if (state==0) state = -1;
//       else if (state==-1) state=0;
//       else return true;

//     } 
//     else //there is a crossing of agent i across the current agent's line
//     { // both agents should be close, or their baseline almost collinear, to allow this happen
//       std::cout<<termcolor::bold<<termcolor::red<<"a rare condition! both crossing happens!"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";     
//       std::cout<<termcolor::bold<<termcolor::blue<<"d1 is"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<d1<<"\n";     
//       std::cout<<termcolor::bold<<termcolor::blue<<"d2 is"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<d2<<"\n";                 
//     }

//   }

//   if (state * c2 * d2 <0) //check anchor point cable length constraint
//   {
//     if ((cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm())<0)
//     {
//       std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";           
//       std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";     
//       std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint!"<<"\n";
//       std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
//       cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm()<<"\n";     
//       return true;      
//     }

//   }
//   return false;
// }

bool eu::checkEntangleUpdate(int& state, int& state2, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik,
                               Eigen::Vector2d& pikplus1, Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength,
                               double maxheight)
{
  // std::cout<<termcolor::bold<<termcolor::green<<"checking entangling condition!"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pk is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pk<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pkplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pkplus1<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pik is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pik<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"state is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<state<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pbi is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pbi<<"\n";

  double c1 = vectorWedge2(pk, pik, pbi);
  double c2 = vectorWedge2(pkplus1, pikplus1, pbi);
  // double d1 = vectorWedge2(pk, pik, pb);
  double d2 = vectorWedge2(pkplus1, pikplus1, pb);
  double e2 = vectorWedge2(pkplus1, pb, pbi);
  double f1 = vectorWedge2(pb, pik, pbi);
  double f2 = vectorWedge2(pb, pikplus1, pbi);
  // bool near_colinear = false;
  // double cond1, cond2;
  if ((f1*f2)<0)
  { 
    // double pikplus1_pkplus1 = (pikplus1-pkplus1).norm();
    // double pb_pkplus1 = (pb-pkplus1).norm();
    // double pbi_pkplus1 = (pbi-pkplus1).norm();

    // if (c1/ ((pik-pk).norm()*(pbi-pk).norm())<0.08 || 
    //     c2/ (pikplus1_pkplus1 * pbi_pkplus1)<0.08||
    //     d2/ (pikplus1_pkplus1 * pb_pkplus1)<0.08 ||
    //     e2/ (pb_pkplus1 * pbi_pkplus1)<0.08)
    // {
    //   near_colinear = true;
    //   cond1 = (pkplus1(0)-pb(0))/(pbi(0)-pb(0))>1.0? 1.0: -1.0;
    //   cond2 = (pkplus1(0)-pikplus1(0))/(pbi(0)-pb(0))>1.0? 1.0: -1.0;
    // }
    // else
    // {
    //   cond1 = e2;
    //   cond2 = d2;
    // }

    if (state ==0) //checking base position of agent i
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and its base!"<<"\n";
      if (e2>0) state = 1;
      else state = -1;
    }
    else if (state * e2 > 0) state = 0;
    // else return true; //else fo nothing for base checking
    if (state2 ==0) //checking position of agent i
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and its base!"<<"\n";
      if (d2>0) state2 = 1;
      else state2 = -1;
    }
    else if (state2 * d2 > 0) state2 = 0;
    // else return true;   //else fo nothing for base checking      


  }  

  if ((c1*c2)<0)
  {

    if (state ==0) //crossing in between agent i and its base
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and its base!"<<"\n";
      if (e2>0) state = 1;
      else state = -1;
    }
    else if (state * e2 > 0) state = 0;
    else
    {
      std::cout<<termcolor::bold<<termcolor::red<<"state * e2 = "<<state * e2<<"\n";
      return true;
    } 
    if (state2 ==0) //crossing in between agent i and its base
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"crossing in between agent i and its base!"<<"\n";
      if (d2>0) state2 = 1;
      else state2 = -1;
    }
    else if (state2 * d2 > 0) state2 = 0;
    else 
    {
      std::cout<<termcolor::bold<<termcolor::red<<"state2 * d2 = "<<state2 * d2<<"\n";
      return true;    
    }

  }

  if (state * e2 <-1e-9)
  {
    if ((cableLength-1.2*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm())<0)
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint for base!"<<"\n";
      std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
      cableLength-1.2*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm()<<"\n";     
      return true;      
    }
    else return false; //when base anchor constraint is active, no need to consider agent anchor constraint    
  }
  if (state2 * d2 <-1e-9 && d2 * e2<0) //check anchor point cable length constraint
  {
    if ((cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm())<0)
    {
      // std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";             
      std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint! for pikplus1"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";         
      std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
      cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm()<<"\n";     
      return true;      
    }

  }
  return false;
}

bool eu::checkEntangleUpdate2(int& state, int& state2, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik,
                               Eigen::Vector2d& pikplus1, Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength,
                               double maxheight)
{
  // std::cout<<termcolor::bold<<termcolor::green<<"checking entangling condition!"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pk is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pk<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pkplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pkplus1<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pik is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pik<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";

  Eigen::Vector2d pik_pb;
  Eigen::Vector2d pbi_pb;
  Eigen::Vector2d pik_pk;
  Eigen::Vector2d pbi_pk;

  double c1 = vectorWedge2(pk, pik, pbi, pik_pk, pbi_pk);
  double c2 = vectorWedge2(pkplus1, pikplus1, pbi);
  double d1 = vectorWedge2(pk, pik, pb);
  // double d2 = vectorWedge2(pkplus1, pikplus1, pb);
  double e1 = vectorWedge2(pk, pbi, pb);

  double f1 = vectorWedge2(pb, pik, pbi, pik_pb, pbi_pb);
  double f2 = vectorWedge2(pb, pikplus1, pbi);

  if ((f1*f2)<0)
  { 
    double a;
    if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
    {
      a = pik_pb(1)/pbi_pb(1);
    }
    else
    {
      a = pik_pb(0)/pbi_pb(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      if (state == 0) 
      {
        if (d1>0) state = 1;
        else state = -1;
        if (e1>0) state2 = 1;
        else state2 = -1;
      }
      else if (state==1 || state ==-1) state =0;
    }
    else if (a<1) //crossing beyond agent i
    {
      if (state == 0) 
      {
        state =2;
        if (e1>0) state2 = 1;
        else state2 = -1;
      }        
      else if (state==2) state =0;      
    }
    else //crossing beyond base i
    {
      if (state == 0)
      {
        state =3;  
        if (e1>0) state2 = 1;
        else state2 = -1;
      } 
      else if (state==3) state =0;            
    }

  }  

  if ((c1*c2)<0)
  {
    double a;
    if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
    {
      a = pik_pk(1)/pbi_pk(1);
    }
    else
    {
      a = pik_pk(0)/pbi_pk(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      if (state == 0) 
      {
        if (d1>0) state = 1;
        else state = -1;
        if (e1>0) state2 = 1;
        else state2 = -1;
      }
      else if (state==1 || state ==-1) state =0;
      else return true;
    }
    else if (a<1) //crossing beyond agent i
    {
      if (state == 0) 
      {
        state =2;
        if (e1>0) state2 = 1;
        else state2 = -1;
      }
      else if (state==2) state =0;   
      else return true;   
    }
    else  //crossing beyong base i
    {
      if (state == 0)
      {
        state =3;  
        if (e1>0) state2 = 1;
        else state2 = -1;
      } 
      else if (state==3) state =0;     
      else return true;       
    }

  }

  if (state != 0 && e1 * state2 < -1e-9) //anchor with base i
  {
    // if ((cableLength-1.2*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm())<0)
    // {
    //   // std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint for base!"<<"\n";
    //   // std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
    //   // cableLength-1.2*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm()<<"\n";     
    //   return true;      
    // }
    // else return false; //when base anchor constraint is active, no need to consider agent anchor constraint    
  }
  
  if ((state == 1 || state == -1) && state * d1 <-1e-9 ) //anchor with the agent i
  {
    // if ((cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm())<0)
    // {
    //   // std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
    //   // std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";             
    //   // std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint! for pikplus1"<<"\n";
    //   // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
    //   // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";         
    //   // std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
    //   // cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm()<<"\n";     
    //   return true;      
    // }
  }

  return false;
}

int eu::entangleStatesUpdate3(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                             const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                             const Eigen::Vector2d& pbi, double& d1, double& e1)
{
  // std::cout<<termcolor::bold<<termcolor::green<<"checking entangling condition!"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pk is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pk<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pkplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pkplus1<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pik is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pik<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
  // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";

  Eigen::Vector2d pik_pb;
  Eigen::Vector2d pbi_pb;
  Eigen::Vector2d pik_pk;
  Eigen::Vector2d pbi_pk;

  double c1 = vectorWedge2(pk, pik, pbi, pik_pk, pbi_pk);
  double c2 = vectorWedge2(pkplus1, pikplus1, pbi);
  d1 = vectorWedge2(pk, pik, pb);
  // double d2 = vectorWedge2(pkplus1, pikplus1, pb);
  e1 = vectorWedge2(pk, pbi, pb);

  double f1 = vectorWedge2(pb, pik, pbi, pik_pb, pbi_pb);
  double f2 = vectorWedge2(pb, pikplus1, pbi);

  if ((f1*f2)<0)
  { 
    double a;
    if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
    {
      a = pik_pb(1)/pbi_pb(1);
    }
    else
    {
      a = pik_pb(0)/pbi_pb(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      if (state == 0) 
      {
        if (d1>0) state = 1;
        else state = -1;
        // state2 = 0;
      }
      else if (state==1 || state ==-1) state =0;
    }
    else if (a<1) //crossing beyond agent i
    {
      if (state == 0) 
      {
        if (d1>0) state =2;
        else state = -2;
        // state2 = 0;
      }        
      else if (state==2 || state ==-2) state =0;      
    }
    else //crossing beyond base i
    {
      if (state == 0)
      {
        if (d1>0) state =3;  
        else state = -3;
        // state2 = 0;
      } 
      else if (state==3 || state ==-3) state =0;            
    }

  }  

  if ((c1*c2)<0)
  {
    double a;
    if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
    {
      a = pik_pk(1)/pbi_pk(1);
    }
    else
    {
      a = pik_pk(0)/pbi_pk(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      if (state3 == 0)
      {
        if (state == 0) 
        {
          if (d1>0) state = 1;
          else state = -1;
          if (e1>0) state2 = 1;
          else state2 = -1;
        }
        else if (state==1 || state ==-1) 
        {
          state =0;
          state2 = 0;
        }
        else
        {
          state3 = 1;
          return RISK_ENTANGLING;
        }         
      }
      else if (state3==1)
      {
        state3 = 0;
      }  
      else   // state3 not 0 or 1
      {
        return ENTANGLED;        
      }    
    }
    else if (a<1) //crossing beyond agent i
    {
      if (state3 == 0)
      {      
        if (state == 0) 
        {
          if (d1>0) state = 2;
          else state = -2;
          if (e1>0) state2 = 1;
          else state2 = -1;
        }
        else if (state==2 || state ==-2) 
        {
          state =0;
          state2 = 0;
        }   
        else 
        {
          state3 = 2;
          return RISK_ENTANGLING;   
        }
      }
      else if (state3==2)
      {
        state3 = 0;
      }
      else return ENTANGLED;
    }
    else  //crossing beyong base i
    {
      if (state3 == 0)
      {         
        if (state == 0)
        {
          if (d1>0) state =3;  
          else state = -3;
          if (e1>0) state2 = 1;
          else state2 = -1;
        } 
        else if (state==3 || state == -3) 
        {
          state =0;
          state2 = 0;
        }     
        else 
        {
          state3 = 3;
          return RISK_ENTANGLING;  
        } 
      }
      else if (state3 == 3)
      {
        state3 = 0;        
      }
      else return ENTANGLED;
    }

  }

  return NOT_ENTANGLED;
}

int eu::checkEntangleUpdate3(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                             const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                             const Eigen::Vector2d& pbi, double cableLength, double maxheight)
{

  double d1, e1;

  int result = entangleStatesUpdate3(state, state2, state3, pk, pkplus1, pik, pikplus1, pb, pbi, d1, e1);

  if (result == RISK_ENTANGLING || result == ENTANGLED) return result;

  if (state != 0 && e1 * state2 < -1e-9) //anchor with base i
  {
    if ((cableLength-1.1*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm())<0)
    {
      // std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint for base!"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
      // cableLength-1.2*maxheight-(pbi-pb).norm()-(pkplus1-pbi).norm()<<"\n";     
      return true;      
    }
    else return false; //when base anchor constraint is active, no need to consider agent anchor constraint    
  }
  
  if ((state == 1 || state == -1) && state * d1 <-1e-9 ) //anchor with the agent i
  {
    if ((cableLength-1.1*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm())<0)
    {
      // std::cout<<termcolor::bold<<termcolor::blue<<"pb is"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<pb<<"\n";             
      // std::cout<<termcolor::bold<<termcolor::red<<"Not satisfying anchor point constraint! for pikplus1"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<"pikplus1 is"<<"\n";
      // std::cout<<termcolor::bold<<termcolor::blue<<pikplus1<<"\n";         
      // std::cout<<termcolor::bold<<termcolor::blue<<"checking value is is"<<
      // cableLength-1.2*maxheight-(pikplus1-pb).norm()-(pkplus1-pikplus1).norm()<<"\n";     
      return true;      
    }
  }

  return NOT_ENTANGLED;
}

int eu::checkEntangleStaticUpdate(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, 
                                  const Eigen::Vector2d& pkplus1, const Eigen::Vector2d& pik, 
                                  const Eigen::Vector2d& pb, const Eigen::Vector2d& pbi, 
                                  double cableLength, double maxheight)
{

  double d1, e1;

  int result = entangleStatesUpdate3(state, state2, state3, pk, pkplus1, pik, pik, pb, pbi, d1, e1);

  if (result == RISK_ENTANGLING || result == ENTANGLED) return result;

  if (state3 == 2)
  {
    if (1.2 * (pik-pbi).norm() + (pik-pkplus1).norm() + (pb-pbi).norm() + maxheight > cableLength)
      return ENTANGLED;
  }
  else if (state3 == 3)
  {
    if (1.2 * (pik-pbi).norm() + (pbi-pkplus1).norm() + (pb-pik).norm() + maxheight > cableLength)
      return ENTANGLED;
  }
  else if (e1 * state2 < -1e-9) 
  {
    double lla = (pb-pbi).norm() + (pkplus1-pbi).norm();
    if (d1 * state < -1e-9)
    {
      double llb = (pb-pik).norm() + (pkplus1-pik).norm();
      if (lla < llb) return checkMaxPointAnchor(llb, cableLength, maxheight);
    }    

    return checkMaxPointAnchor(lla, cableLength, maxheight);
  }
  else if ( state * d1 <-1e-9 ) //anchor with the agent i
  {
    double llb = (pb-pik).norm() + (pkplus1-pik).norm();
    return checkMaxPointAnchor(llb, cableLength, maxheight);
  }

  return NOT_ENTANGLED;
}

int eu::checkMaxPointAnchor(double ll, double cableLength, double maxheight)
{
  if ((cableLength-1.1*maxheight-ll)<0)
  {   
    return ENTANGLED;      
  }  
  else return NOT_ENTANGLED;
}

void eu::entangleHSigUpdate3(int id, std::vector<Eigen::Vector2i>& alphas, std::vector<Eigen::Vector2d>& betas, 
                            int& active_case_num, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                             const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                             const Eigen::Vector2d& pbi)
{
  int active_case_num_old = active_case_num;

  Eigen::Vector2d pik_pb;
  Eigen::Vector2d pbi_pb;
  Eigen::Vector2d pik_pk;
  Eigen::Vector2d pbi_pk;

  double c1 = vectorWedge2(pk, pik, pbi, pik_pk, pbi_pk);
  double c2 = vectorWedge2(pkplus1, pikplus1, pbi);
  double d1 = vectorWedge2(pk, pik, pb);
  // double d2 = vectorWedge2(pkplus1, pikplus1, pb);
  double e1 = vectorWedge2(pk, pbi, pb);

  double f1 = vectorWedge2(pb, pik, pbi, pik_pb, pbi_pb);
  double f2 = vectorWedge2(pb, pikplus1, pbi);

  if ((f1*f2)<0)
  { 
    double a;
    if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
    {
      a = pik_pb(1)/pbi_pb(1);
    }
    else
    {
      a = pik_pb(0)/pbi_pb(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
    }
    else if (a<1) //crossing beyond agent i
    {
      addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
    }
    else //crossing beyond base i
    {
      addOrRemoveLastElement(id, 3, alphas, betas, active_case_num, d1, e1);
    }

  }  

  if ((c1*c2)<0)
  {
    double a;
    if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
    {
      a = pik_pk(1)/pbi_pk(1);
    }
    else
    {
      a = pik_pk(0)/pbi_pk(0);
    }

    if (a<0) //crossing in between agent i and its base
    {
      addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
    }
    else if (a<1) //crossing beyond agent i
    {
      addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
    }
    else  //crossing beyong base i
    {
      addOrRemoveLastElement(id, 3, alphas, betas, active_case_num, d1, e1);
    }

  }

  // if (active_case_num<=1) return NOT_ENTANGLED;
  // else if (active_case_num==2) return RISK_ENTANGLING;
  // else return ENTANGLED;
  // return NOT_ENTANGLED;
}

void eu::addOrRemoveLastElement(int agent_id, int case_id, std::vector<Eigen::Vector2i>& alphas, 
                                std::vector<Eigen::Vector2d>& betas, int& active_case_num, double d1, double e1)
{
    if (!alphas.empty())
    {
      for (int i=alphas.size()-1; i>=0; i--)
      {
        if (alphas[i](0)== agent_id && alphas[i](1) == case_id)
        {
          // alphas.pop_back();
          // betas.pop_back();
          alphas.erase (alphas.begin()+i);
          betas.erase (betas.begin()+i);
          active_case_num -= 1;        
          return;  
        }

        // if (alphas[i](0)<=4 && alphas[i](1) ==1 || alphas[i](0)== agent_id) break; //as long as there is one non-static obs
        if (agent_id<=4 && case_id!=1)
        {
          if (alphas[i](0)<=4 && alphas[i](1) ==1 || alphas[i](0)== agent_id) break;
        }
        else if (agent_id>4)
        {
          if (alphas[i](0)<=4 && alphas[i](1) ==1 || alphas[i](0) >4) break;
        }
        else break;
      }

    }

    alphas.push_back(Eigen::Vector2i(agent_id, case_id));
    betas.push_back(Eigen::Vector2d(d1, e1));
    active_case_num += 1;
}

void eu::shrinkAlphaBetas(std::vector<Eigen::Vector2i>& alphas, std::vector<Eigen::Vector2d>& betas, 
                               std::vector<int>& active_cases, int num_of_agents)
{
  for (int i=0; i<active_cases.size(); i++)
  {
    if (active_cases[i]>1)
    {
      for (int j=0; j<alphas.size(); j++) // check the entry 
      {
        if (alphas[j](0) == i+1)
        {
          // int index_to_check = j;
          for (int k = j+1; k<alphas.size(); k++)
          {
            if (alphas[j] == alphas[k])
            {
              alphas.erase (alphas.begin()+k); //remove k first
              betas.erase (betas.begin()+k);
              alphas.erase (alphas.begin()+j);
              betas.erase (betas.begin()+j);
              active_cases[i] -= 2;      

              // if a removal happens, redo using the new order
              shrinkAlphaBetas(alphas, betas, active_cases, num_of_agents); 
              return;
            }
            if (alphas[k](0) < num_of_agents && alphas[k](1) == 1) break; //is a dynamic agent 
            if (alphas[k](0) >= num_of_agents && alphas[j](0) >= num_of_agents) break;            
          }
        }
      }
    }
  }
}

void eu::entangleHSigToAddAgentInd(std::vector<Eigen::Vector2i>& alphasToAdd, const Eigen::Vector2d& pk, 
                                   const Eigen::Vector2d& pkplus1, const Eigen::Vector2d& pik, 
                                   Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                                   std::vector<Eigen::Vector2d>& bendpts, std::vector<Eigen::Vector2d>& bendpts_prev, 
                                   int agent_id)
{
  if (bendpts_prev.size()==bendpts.size())
  {
    entangleHSigToAddAgentInd(alphasToAdd, pk, pkplus1, pik, pikplus1, pb, bendpts, agent_id);
    return;
  }

  bool have_base_addition = false;
  if (bendpts.size() ==0 || bendpts_prev.size()==0) //separate the one segement case from >1 seg case
  {
    std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: bendpts empty!"<<"\n";
  }
  else if (bendpts.size() < bendpts_prev.size()) //getting released from bendpoint
  {
    for (int i=0; i< bendpts.size(); i++)
    { 
      Eigen::Vector2d pik_pk, pik_pk_prime;
      Eigen::Vector2d pbi_pk, pbi_pk_prime;
      double c1;
      double c2;
      double c1_prime;  

      if (i != bendpts.size()-1) 
      {
        c1 = vectorWedge2(pk, bendpts[i+1], bendpts[i], pik_pk, pbi_pk);
        c2 = vectorWedge2(pkplus1, bendpts[i+1], bendpts[i]);  
      }
      else
      { 
        c1 = vectorWedge2(pk, pik, bendpts_prev.back(), pik_pk, pbi_pk);
        c2 = vectorWedge2(pkplus1, pikplus1, bendpts[i]);   
        c1_prime = vectorWedge2(pk, bendpts_prev.back(), bendpts_prev.end()[-2], 
                                pik_pk_prime, pbi_pk_prime);        
      }

      if (i == bendpts.size()-1)  //only for the last seg, check for base update
      {
        Eigen::Vector2d pik_pb, pbi_pb;      
        double f1 = vectorWedge2(pb, pik, bendpts_prev.back());
        double f2 = vectorWedge2(pb, pikplus1, bendpts[i], pik_pb, pbi_pb);        
        double f1_prime; 
        if (i==0) f1_prime = vectorWedge2(pb, bendpts_prev.back(), bendpts_prev.end()[-2]);

        if ((f1*f2)<0)  
        { 
          double a;
          if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
          {
            a = pik_pb(1)/pbi_pb(1);
          }
          else
          {
            a = pik_pb(0)/pbi_pb(0);
          }

          if (a<0) //crossing in between agent i and its base
          {
            // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
          }
          else if (a<1) //crossing beyond agent i
          {
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
            // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
          }
          else if (i==0)//crossing beyond base i
          {
            //only if last seg is also the first seg, i.e. i=0
            // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));        
          }
          have_base_addition = true;
        }
        if (i==0 && (f1_prime*f2)<0) //only for the case of only one seg, check for this
        {
          double a;
          if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
          {
            a = pik_pb(1)/pbi_pb(1);
          }
          else
          {
            a = pik_pb(0)/pbi_pb(0);
          }

          if (a<0) //crossing in between agent i and its base
          {
            // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
          }
          else if (a<1) //not to be checked in this case
          {
            // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
            // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
          }
          else//crossing beyond base i
          {
            //only if last seg is also the first seg, i.e. i=0
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));     
            std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: stop4!"<<"\n";
            std::cout<<termcolor::bold<<termcolor::red<<"i is "<<i<<" bendptsize is "<<bendpts.size()<<"\n";

            exit(-1);
   
          }
          have_base_addition = true;          
        }                  
      }

      bool added_a_inbtw_crossing = false;
      if ((c1*c2)<0)
      {
        double a;
        if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
        {
          a = pik_pk(1)/pbi_pk(1);
        }
        else
        {
          a = pik_pk(0)/pbi_pk(0);
        }

        if (a<0) //crossing in between agent i and its base
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, i+2));
          added_a_inbtw_crossing = true;
        }
        else if (a<1 && i == bendpts.size()-1) //crossing beyond agent i, only when checking last seg
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
        }
        else if (a>=1 && i == 0) //crossing beyond base i, not in this case because it is check in the next one
        {
          // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));
        }
      }
      if (i == bendpts.size()-1 && c1_prime*c2<0)
      {
        double a;
        if (abs(pik_pk_prime(1)*pbi_pk_prime(1))>abs(pik_pk_prime(0)*pbi_pk_prime(0)))
        {
          a = pik_pk_prime(1)/pbi_pk_prime(1);
        }
        else
        {
          a = pik_pk_prime(0)/pbi_pk_prime(0);
        }

        if (a<0 && !added_a_inbtw_crossing) //crossing in between agent i and its base
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, i+2));
        }
        else if (a<1 && i == bendpts.size()-1) //crossing beyond agent i, not possible in this case
        {
          // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
        }
        else if (a>=1 && i == 0) //crossing beyond base i, only possible when there is one seg
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));
            std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: stop3!"<<"\n";
            std::cout<<termcolor::bold<<termcolor::red<<"i is "<<i<<" bendptsize is "<<bendpts.size()<<"\n";

            exit(-1);          
        }        
      }
    }    
  }
  else // getting added a bendpoint
  {
    for (int i=0; i< bendpts.size(); i++)
    { 
      Eigen::Vector2d pik_pk;
      Eigen::Vector2d pbi_pk;
      double c1;
      double c2;

      if (i == bendpts.size()-1) 
      {
        c1 = vectorWedge2(pk, pik, bendpts_prev.back());
        c2 = vectorWedge2(pkplus1, pikplus1, bendpts[i], pik_pk, pbi_pk);  
      }
      else if (i == bendpts.size()-2)
      {
        c1 = vectorWedge2(pk, pik, bendpts_prev.back());
        c2 = vectorWedge2(pkplus1, bendpts[i+1], bendpts[i], pik_pk, pbi_pk);          
      }
      else
      {  
        c1 = vectorWedge2(pk, bendpts[i+1], bendpts[i], pik_pk, pbi_pk);
        c2 = vectorWedge2(pkplus1, bendpts[i+1], bendpts[i]);       
      }

      if (i == bendpts.size()-1)  //only for the last seg, check for base update
      {
        Eigen::Vector2d pik_pb, pbi_pb;      
        double f1 = vectorWedge2(pb, pik, bendpts_prev.back());
        double f2 = vectorWedge2(pb, pikplus1, bendpts[i], pik_pb, pbi_pb);        

        if ((f1*f2)<0)  
        { 
          double a;
          if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
          {
            a = pik_pb(1)/pbi_pb(1);
          }
          else
          {
            a = pik_pb(0)/pbi_pb(0);
          }

          if (a<0) //crossing in between agent i and its base
          {
            // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
          }
          else if (a<1) //crossing beyond agent i
          {
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
            // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
          }
          else if (i==0)//crossing beyond base i
          {
            //only if last seg is also the first seg, i.e. i=0
            // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));        
          }
          have_base_addition = true;
        }               
      }
      if (i == 0 && bendpts.size()==2)  //only when exactly two segments, which means first segment is moved
      {     
        Eigen::Vector2d pik_pb, pbi_pb;      
        double f1 = vectorWedge2(pb, pik, bendpts_prev.back());
        double f2 = vectorWedge2(pb, bendpts[i+1], bendpts[i], pik_pb, pbi_pb);  

        if (f1*f2<0) //only for the case of only one seg, check for this
        {
          double a;
          if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
          {
            a = pik_pb(1)/pbi_pb(1);
          }
          else
          {
            a = pik_pb(0)/pbi_pb(0);
          }

          if (a<0) //crossing in between agent i and its base
          {
            // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
          }
          else if (a<1) //not to be checked in this case
          {
            // alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
            // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
          }
          else//crossing beyond base i
          {
            //only if last seg is also the first seg, i.e. i=0
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0)); 
            std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: stop2!"<<"\n";
            std::cout<<termcolor::bold<<termcolor::red<<"i is "<<i<<" bendptsize is "<<bendpts.size()<<"\n";

            exit(-1);                   
          }
          have_base_addition = true;          
        }
      }   
      if ((c1*c2)<0)
      {
        double a;
        if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
        {
          a = pik_pk(1)/pbi_pk(1);
        }
        else
        {
          a = pik_pk(0)/pbi_pk(0);
        }

        if (a<0) //crossing in between agent i and its base
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, i+2));
        }
        else if (a<1 && i == bendpts.size()-1) //crossing beyond agent i, only when checking last seg
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
        }
        else if (a>=1 && i == 0) //crossing beyond base i
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));
            std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: stop1!"<<"\n";
            std::cout<<termcolor::bold<<termcolor::red<<"i is "<<i<<" bendptsize is "<<bendpts.size()<<"\n";

            exit(-1);          
        }
      }
    } 
  }  

  if (have_base_addition)
  {
    if (alphasToAdd.size()>=2 && alphasToAdd.back() == alphasToAdd.end()[-2])
    {
      alphasToAdd.erase(alphasToAdd.end()-2, alphasToAdd.end());
    }
  }
}

void eu::entangleHSigToAddAgentInd(std::vector<Eigen::Vector2i>& alphasToAdd,  
                            const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                            const Eigen::Vector2d& pik, Eigen::Vector2d& pikplus1,
                            const Eigen::Vector2d& pb, std::vector<Eigen::Vector2d>& bendpts, int agent_id)
{
  bool have_base_addition = false;

  if (bendpts.size() >= 1)
  {
    for (int i=0; i< bendpts.size(); i++)
    { 
      Eigen::Vector2d pik_pk;
      Eigen::Vector2d pbi_pk;
      double c1;
      double c2;  

      if (i != bendpts.size()-1) 
      {
        c1 = vectorWedge2(pk, bendpts[i+1], bendpts[i], pik_pk, pbi_pk);
        c2 = vectorWedge2(pkplus1, bendpts[i+1], bendpts[i]);  
      }
      else
      { 
        c1 = vectorWedge2(pk, pik, bendpts[i], pik_pk, pbi_pk);
        c2 = vectorWedge2(pkplus1, pikplus1, bendpts[i]);  
      }

      if (i == bendpts.size()-1)  //only for the last seg, check for base update
      {
        Eigen::Vector2d pik_pb, pbi_pb;
        double f1 = vectorWedge2(pb, pik, bendpts[i], pik_pb, pbi_pb);
        double f2 = vectorWedge2(pb, pikplus1, bendpts[i]);        

        if ((f1*f2)<0)
        { 
          double a;
          if (abs(pik_pb(1)*pbi_pb(1))>abs(pik_pb(0)*pbi_pb(0)))
          {
            a = pik_pb(1)/pbi_pb(1);
          }
          else
          {
            a = pik_pb(0)/pbi_pb(0);
          }

          if (a<0) //crossing in between agent i and its base
          {
            // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
          }
          else if (a<1) //crossing beyond agent i
          {
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
            // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
          }
          else if (i==0)//crossing beyond base i
          {
            //only if last seg is also the first seg, i.e. i=0
            alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));        
          }
          have_base_addition = true;
        }          
      }

      if ((c1*c2)<0)
      {
        double a;
        if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
        {
          a = pik_pk(1)/pbi_pk(1);
        }
        else
        {
          a = pik_pk(0)/pbi_pk(0);
        }

        if (a<0) //crossing in between agent i and its base
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, i+2));
        }
        else if (a<1 && i == bendpts.size()-1) //crossing beyond agent i, only when checking last seg
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 1));
        }
        else if (a>=1 && i == 0) //crossing beyond base i, only when checking first seg
        {
          alphasToAdd.push_back(Eigen::Vector2i(agent_id, 0));
        }
      }
    }    
  }
  else std::cout<<termcolor::bold<<termcolor::red<<"entangleHSigToAddAgentInd: bendpts empty!"<<"\n";

  if (have_base_addition)
  {
    if (alphasToAdd.size()>=2 && alphasToAdd.back() == alphasToAdd.end()[-2])
    {
      alphasToAdd.erase(alphasToAdd.end()-2, alphasToAdd.end());
    }
  }  
}


void eu::entangleHSigToAddStatic(std::vector<Eigen::Vector2i>& alphasToAdd,  
                            const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                            std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents)
{
  for (int i=0; i<staticObsRep.size(); i++)
  {
    Eigen::Vector2d pik_pk;
    Eigen::Vector2d pbi_pk;

    Eigen::Vector2d pik = staticObsRep[i].col(1);
    Eigen::Vector2d pbi = staticObsRep[i].col(0);

    double c1 = vectorWedge2(pk, pik, pbi, pik_pk, pbi_pk);
    double c2 = vectorWedge2(pkplus1, pik, pbi);    

    if ((c1*c2)<0)
    {
      double a;
      if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
      {
        a = pik_pk(1)/pbi_pk(1);
      }
      else
      {
        a = pik_pk(0)/pbi_pk(0);
      }

      if (a<0) //this should not happen for static obst
      {
        // addOrRemoveLastElement(id, 1, alphas, betas, active_case_num, d1, e1);
        // alphasToAdd.push_back(Eigen::Vector2i(num_of_agents+i+1, 1));
      }
      else if (a<1) //crossing beyond agent i
      {
        // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
        alphasToAdd.push_back(Eigen::Vector2i(num_of_agents+i+1, 1));
      }
      else  //crossing beyong base i
      {
        // addOrRemoveLastElement(id, 3, alphas, betas, active_case_num, d1, e1);
        alphasToAdd.push_back(Eigen::Vector2i(num_of_agents+i+1, 0));
      }

    }

  }
}

void eu::updateHSig_orig(ent_state_orig& entangle_state,  
                        const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                        std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep)
{
  std::vector<Eigen::Vector2i> alphasToAdd;
  for (int i=0; i<staticObsRep.size(); i++)
  {
    Eigen::Vector2d pik_pk;
    Eigen::Vector2d pbi_pk;

    Eigen::Vector2d pik = staticObsRep[i].col(1);
    Eigen::Vector2d pbi = staticObsRep[i].col(0);

    double c1 = vectorWedge2(pk, pik, pbi, pik_pk, pbi_pk);
    double c2 = vectorWedge2(pkplus1, pik, pbi);    

    if ((c1*c2)<0)
    {
      double a;
      if (abs(pik_pk(1)*pbi_pk(1))>abs(pik_pk(0)*pbi_pk(0)))
      {
        a = pik_pk(1)/pbi_pk(1);
      }
      else
      {
        a = pik_pk(0)/pbi_pk(0);
      }

      if (a<0) //this should not happen for static obst
      {

      }
      else if (a<1) //crossing beyond agent i
      {
        // addOrRemoveLastElement(id, 2, alphas, betas, active_case_num, d1, e1);
        if (c1<0) alphasToAdd.push_back(Eigen::Vector2i(i+1, 1));
        else alphasToAdd.push_back(Eigen::Vector2i(i+1, -1));

      }
      else  //crossing beyong base i
      {

      }
    }
  }

  bool have_cancellation = true;
  while (have_cancellation)
  {
    for (int i=0; i<alphasToAdd.size(); i++)
    {
      for (int j=entangle_state.hsig.size()-1; j>=0; j--) //should only check until last bend point
      {
        if (entangle_state.hsig[j](0) == alphasToAdd[i](0) &&
            entangle_state.hsig[j](1) != alphasToAdd[i](1)) 
        {
          alphasToAdd.erase (alphasToAdd.begin()+i);
          entangle_state.hsig.erase (entangle_state.hsig.begin()+j);
          
          have_cancellation = true;  
          goto endwhile;
        }
        break;        
      } 
    }
    have_cancellation = false;
    endwhile:;
  }

  if (alphasToAdd.empty()) return; //nothing to add

  for (int i=0; i<alphasToAdd.size(); i++)
  {
    entangle_state.hsig.push_back(alphasToAdd[i]);
  }    
}

void eu::updateContPts_orig(std::vector<Eigen::Vector2d>& contPts,  
                            std::vector<Eigen::Matrix<double, 2, -1>>& convexHullOfStaticObs)
{
  double interval = 0.2;
  std::vector<Eigen::Vector2d> newContPts;
  newContPts.push_back(contPts[0]);
  Eigen::Vector2d currentContPt = contPts[0];
  Eigen::Vector2d pt_prev = currentContPt;
  for (int i=0; i<contPts.size()-1; i++)
  {
    int ratio = static_cast <int> (std::floor((contPts[i]-contPts[i+1]).norm()/interval));
    for (int j=1; j<=ratio; j++)
    {
      Eigen::Vector2d pt = contPts[i] + (double)j* interval*
                           (contPts[i+1]-contPts[i])/(contPts[i]-contPts[i+1]).norm();
      Eigen::Matrix<double, 2, 2> line;
      line << currentContPt, pt;

      for (int k=0; k<convexHullOfStaticObs.size(); k++)
      {
        if (gjk::collision(convexHullOfStaticObs[k], line))
        {
          currentContPt = pt_prev;
          newContPts.push_back(currentContPt);
          break;          
        }
      }
      pt_prev = pt;
    }

    Eigen::Matrix<double, 2, 2> line;
    line << currentContPt, contPts[i+1];
    for (int k=0; k<convexHullOfStaticObs.size(); k++)
    {
      if (gjk::collision(convexHullOfStaticObs[k], line))
      {
        currentContPt = pt_prev;
        newContPts.push_back(currentContPt);
        break;          
      }
    }
    pt_prev = contPts[i+1];    
  }
  newContPts.push_back(contPts.back());
  contPts = newContPts;
}
void eu::addAlphaBetaToList(std::vector<Eigen::Vector2i>& alphasToAdd, ent_state& entangle_state, 
                            const Eigen::Vector2d& pk, std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& pb_self,
                            std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents,
                            std::vector<std::vector<Eigen::Vector2d>>& bendpts)
{
  int have_cancellation = true;
  //firstly, try to (sequentially) cancel the bend point if possible
  // while (have_cancellation && !entangle_state.bendPointsIdx.empty())
  // {
  //   for (int i=0; i<alphasToAdd.size(); i++)
  //   {
  //     if (alphasToAdd[i] != entangle_state.alphas[entangle_state.bendPointsIdx.back()]) 
  //       continue; //check if it can potentially cancel the bendpoint

  //     for (int j=entangle_state.alphas.size()-1; j>=0; j--) //check from the last point in list
  //     {
  //       if (entangle_state.alphas[j] == alphasToAdd[i]) 
  //       {
  //         if (j != entangle_state.bendPointsIdx.back()) break; //only cancel if it matches the bendpoint
  //         entangle_state.active_cases[alphasToAdd[i](0)-1] -= 1;  
  //         alphasToAdd.erase (alphasToAdd.begin()+i);
  //         entangle_state.alphas.erase (entangle_state.alphas.begin()+j);
  //         entangle_state.betas.erase (entangle_state.betas.begin()+j);

  //         entangle_state.bendPointsIdx.pop_back();
  //         Eigen::Vector2d bp;
  //         getBendPt2d(bp, entangle_state, pb, pb_self, staticObsRep, num_of_agents);          
  //         for (int k=j; j<entangle_state.alphas.size(); j++) //start from the new j (old j has been erased)
  //         { //update the betas for the new bendpoint
  //           entangle_state.betas[j] = calculateBetaForCase(entangle_state.alphas[j], pk, pb, 
  //                                                          bp, staticObsRep, num_of_agents);
  //         }

  //         have_cancellation = true;  
  //         goto endwhile1;
  //       }

  //       if (breakcondition(alphasToAdd[i], entangle_state.alphas[j], num_of_agents, j, 
  //                        entangle_state.bendPointsIdx.back())) break;
  //     }      
  //   }
  //   have_cancellation = false;
  //   endwhile1:;
  // }

  // have_cancellation = true;
  // int _b = entangle_state.bendPointsIdx.empty()? 0 : entangle_state.bendPointsIdx.back();

  //then check cancellation of the non-bend points
  while (have_cancellation)
  {
    int _b = entangle_state.bendPointsIdx.empty()? -1 : entangle_state.bendPointsIdx.back();
    for (int i=0; i<alphasToAdd.size(); i++)
    {
      for (int j=entangle_state.alphas.size()-1; j>=0; j--) //should only check until last bend point
      {
        if (entangle_state.alphas[j] == alphasToAdd[i] ||
          alphasToAdd[i](0)<=num_of_agents && entangle_state.alphas[j](0)==alphasToAdd[i](0) && //same agent 
          alphasToAdd[i](1)>=bendpts[alphasToAdd[i](0)-1].size()+1 && //to be added is last segment of agent
          alphasToAdd[i](1) <entangle_state.alphas[j](1) || //the segment number in list is even larger
          //second last segment but used to be the last seg
          // alphasToAdd[i](0)<=num_of_agents && entangle_state.alphas[j](0)==alphasToAdd[i](0) &&
          // entangle_state.alphas[j](1)>=2 && entangle_state.alphas[j](1)+1==alphasToAdd[i](1) &&
          // alphasToAdd[i](1)==bendpts[alphasToAdd[i](0)-1].size()+1 && j > _b ||
          alphasToAdd[i](0)<=num_of_agents && entangle_state.alphas[j](0)==alphasToAdd[i](0) && //same agent 
          entangle_state.alphas[j](1)>=2 && alphasToAdd[i](1) >=2 &&
          abs(alphasToAdd[i](1) - entangle_state.alphas[j](1)) == 1 && j > _b ) 
        {
          // if (alphasToAdd[i](0)<=num_of_agents && alphasToAdd[i](1)<2)
          // { //if it is agent crossing beyond case, check if there is another same agent crossing in-btw case
          //   for (int k=j-1; k>=0; k--) 
          //   {
          //     if (entangle_state.alphas[k](0) == alphasToAdd[i](0) && 
          //         entangle_state.alphas[k](1)>=2)
          //       goto nextitem; //try to cancel the in-btw case first before beyond case
          //   }
          // }

          entangle_state.active_cases[alphasToAdd[i](0)-1] -= 1;  
          alphasToAdd.erase (alphasToAdd.begin()+i);
          entangle_state.alphas.erase (entangle_state.alphas.begin()+j);
          entangle_state.betas.erase (entangle_state.betas.begin()+j);
          
          if (j == _b)
          {
            entangle_state.bendPointsIdx.pop_back();
            Eigen::Vector2d bp;
            getBendPt2d(bp, entangle_state, pb, pb_self, staticObsRep, num_of_agents);          
            for (int k=j; k<entangle_state.alphas.size(); k++) //start from the new j (old j has been erased)
            { //update the betas for the new bendpoint
              entangle_state.betas[k] = calculateBetaForCase(entangle_state.alphas[k], pk, pb, 
                                                             bp, staticObsRep, num_of_agents);
            }
          }
          else if (j < _b)
          {
            entangle_state.bendPointsIdx.back() = _b - 1;
            for (int k = entangle_state.bendPointsIdx.size()-2; k >=0; k--)
            {
              if (entangle_state.bendPointsIdx[k] >j)
              {
                entangle_state.bendPointsIdx[k] -= 1;
              }
              else break;
            }
          }
          have_cancellation = true;  
          goto endwhile;
        }

        if (breakcondition(alphasToAdd[i], entangle_state.alphas[j], num_of_agents, j, _b, bendpts))
          break;        
      } 
      nextitem:;     
    }
    have_cancellation = false;
    endwhile:;
  }

  if (alphasToAdd.empty()) return; //nothing to add

  Eigen::Vector2d _pb; 
  getBendPt2d(_pb, entangle_state, pb, pb_self, staticObsRep, num_of_agents);

  for (int i=0; i<alphasToAdd.size(); i++)
  {
    entangle_state.alphas.push_back(alphasToAdd[i]);
    entangle_state.active_cases[alphasToAdd[i](0)-1] += 1;  

    entangle_state.betas.push_back(
      calculateBetaForCase(alphasToAdd[i], pk, pb, _pb, staticObsRep, num_of_agents));
  }
}

void eu::updateBendPts(ent_state& entangle_state, 
                      const Eigen::Vector2d& pkplus1, std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& pb_self,
                      std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents)
{
  Eigen::Vector2d bp; 
  getBendPt2d(bp, entangle_state, pb, pb_self, staticObsRep, num_of_agents);
  // if there is no bendpt fall back, check if new bendpoint is added
  int idx_new_bend = -1;
  int idx_start = entangle_state.bendPointsIdx.empty()? -1 : entangle_state.bendPointsIdx.back();
  for (int i=idx_start+1; i<entangle_state.alphas.size(); i++)
  {
    double beta = calculateBetaForCase(entangle_state.alphas[i], pkplus1, pb, 
                                       bp, staticObsRep, num_of_agents);
    if (beta * entangle_state.betas[i] <-1e-7)
    {
      idx_new_bend = i;
    }
  }
  if (idx_new_bend>-1) //if new bendpoint is added
  {
    entangle_state.bendPointsIdx.push_back(idx_new_bend);

    Eigen::Vector2d _bp;
    Eigen::Vector2i bend_id = entangle_state.alphas[idx_new_bend];
    if (bend_id(0) <= num_of_agents) _bp = pb[bend_id(0)-1]; //only base of mobile agent possible
    else _bp = staticObsRep[bend_id(0)-num_of_agents-1].col(bend_id(1)); // static obst

    //update the following betas wrt the new bend point
    for (int i=idx_new_bend+1; i<entangle_state.alphas.size(); i++)
    {
      entangle_state.betas[i] = calculateBetaForCase(entangle_state.alphas[i], pkplus1, pb, 
                                       _bp, staticObsRep, num_of_agents);
    }
    return; //if new bendpoint is added, no need to check the old bend points anymore
  }

  // now, check if some existing bend points are released
  // bool have_fallback = false;
  while (!entangle_state.bendPointsIdx.empty())
  { 
    Eigen::Vector2d bp_prev; 
    // get bp_prev
    if (entangle_state.bendPointsIdx.size()==1) bp_prev = pb_self; 
    else
    {
      Eigen::Vector2i bend_id = entangle_state.alphas[    //get second last element
                                entangle_state.bendPointsIdx[entangle_state.bendPointsIdx.size()-2]];

      if (bend_id(0) <= num_of_agents) bp_prev = pb[bend_id(0)-1]; //only base of mobile agent possible
      else bp_prev = staticObsRep[bend_id(0)-num_of_agents-1].col(bend_id(1)); // static obst
    }
    double beta = calculateBetaForCase(entangle_state.alphas[entangle_state.bendPointsIdx.back()], 
                                       pkplus1, pb, bp_prev, staticObsRep, num_of_agents);
    if (beta * entangle_state.betas[entangle_state.bendPointsIdx.back()] > 1e-7)
    {
      // have_fallback = true;
            
      for (int k=entangle_state.bendPointsIdx.back()+1; k<entangle_state.alphas.size(); k++) 
      { //update the betas for the new bendpoint
        entangle_state.betas[k] = calculateBetaForCase(entangle_state.alphas[k], pkplus1, pb, 
                                                       bp_prev, staticObsRep, num_of_agents);
      }
      entangle_state.bendPointsIdx.pop_back();      
    }
    else break;
  }
  // if have fall back, simply return
  // if (have_fallback) return;
}

//return true:  check will stop here, cannot check the next one
//return false: checker can skip current one and continue to check if next one can be cancelled
bool eu::breakcondition(Eigen::Vector2i& alphaToAdd, Eigen::Vector2i& alphaInList, int num_of_agents,
                        int idx_to_check, int idx_last_bend, 
                        std::vector<std::vector<Eigen::Vector2d>>& bendpts) 
{
  // if (alphas[i](0)<=4 && alphas[i](1) ==1 || alphas[i](0)== agent_id) break; //as long as there is one non-static obs
  if (alphaToAdd(0) <= num_of_agents && alphaToAdd(1) >=2) // mobile agents and in-btw case, should not skip 
  {
    if (alphaToAdd(1) <= bendpts[alphaToAdd(0)-1].size()) //not the last seg
    {
      // if (alphaInList(0)<=num_of_agents && alphaInList(1) >=2|| //not skip if in list is agent in-btw crossing case
      //     alphaInList(0) >num_of_agents) //not skpi if in list is static obst
      //  return true; 
      // if (alphaInList(0) == alphaToAdd(0)) return true; // not skip if in list is the same agent but diff case
    }
    else //agent last seg can skip through static obst because it moves a lot
    {
      // if (alphaInList(0)<=num_of_agents && alphaInList(1) >=2) //not skip if in list is agent in-btw crossing case
      //  return true; 
    }

    if (idx_to_check <= idx_last_bend) return true; // not skip if exceeding the current bend list) )
  }
  else if (alphaToAdd(0) <= num_of_agents && alphaToAdd(1) <2) // mobile agents and beyond case, usually can skip
  {
    // if (alphaInList(0) == alphaToAdd(0)) return true; // unless in list is the same agent but diff case
  }
  else if (alphaToAdd(0) > num_of_agents)  //static obst
  {
    if (alphaInList(0) >num_of_agents || // not skip if in list is also static
        idx_to_check <= idx_last_bend) // not skip if exceeding the current bend list) 
        return true; // not skip
    // if (alphaInList(0)<=num_of_agents && alphaInList(1) >= 2 && // inlist is mobile agents and in-btw case AND
    //     alphaInList(1)<=bendpts[alphaInList(0)-1].size()) //inlist is not the last segment crossing
    //   return true;
    //if (alphaInList(0) == alphaToAdd(0)) return true; // unless in list is the same static obst but diff case

  }

  return false;
}

void eu::getBendPt2d(Eigen::Vector2d& bp, eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                     Eigen::Vector2d& pb_self, std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                     int num_of_agents)
{
  if (entangle_state.bendPointsIdx.empty()) bp = pb_self; // bendPointsIdx empty means no intermediate bend point
  else
  {
    Eigen::Vector2i bend_id = entangle_state.alphas[entangle_state.bendPointsIdx.back()];

    if (bend_id(0) > num_of_agents + staticObsRep.size() ||
        bend_id(0)<0 ||
        bend_id(0) > num_of_agents && (bend_id(1)>1 || bend_id(1)<0))
    {
      std::cout<<"something seriously wrong man!\n";
      std::cout<<"entangle_state.betas size: "<<entangle_state.betas.size()<<"\n";
      std::cout<<"entangle_state.alphas size: "<<entangle_state.alphas.size()<<"\n";
      std::cout<<"entangle_state.bendPointsIdx back: "<<entangle_state.bendPointsIdx.back()<<"\n";
      std::cout<<"bend_id: "<<bend_id<<"\n";

    }

    if (bend_id(0) <= num_of_agents && bend_id(0) >=1) 
    {
      bp = pb[bend_id(0)-1]; //only base of mobile agent possible
    }
    else if (bend_id(0) > num_of_agents) 
    {
      bp = staticObsRep[bend_id(0)-num_of_agents-1].col(bend_id(1)); // static obst
    }
  }
}

void eu::getBendPt2dwIdx(Eigen::Vector2d& bp, eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                         std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                         std::vector<Eigen::Vector2d>& staticObsLongestDist,
                         int num_of_agents, int idxFromBeg, double& compensation)
{
  // if (entangle_state.bendPointsIdx.empty()) // bendPointsIdx empty means no intermediate bend point
  // {
  //   std::cout<<termcolor::bold<<termcolor::red<<
  //             "getBendPt2dwIdx: empty entangle_state.bendPointsIdx!"<<"\n";
  //   return;
  // }
  // else
  // {
    Eigen::Vector2i bend_id = entangle_state.alphas[idxFromBeg];

    if (bend_id(0) <= num_of_agents) 
    {
      bp = pb[bend_id(0)-1]; //only base of mobile agent possible
      compensation = 0.0;
    }
    else 
    {
      bp = staticObsRep[bend_id(0)-num_of_agents-1].col(bend_id(1)); // static obst
      compensation = staticObsLongestDist[bend_id(0)-num_of_agents-1](bend_id(1));
    }
  // }
}

double eu::calculateBetaForCase(const Eigen::Vector2i& alphaToAddOrUpdate, const Eigen::Vector2d& pk, 
                            std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& bp, 
                            std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents)
{
  if (alphaToAddOrUpdate(0) <= num_of_agents)
  {
    if (alphaToAddOrUpdate(1)==2 || alphaToAddOrUpdate(1)==0) 
    {
      return 0.0;//vectorWedge2(pk, pb[alphaToAddOrUpdate(0)-1], bp);
    }
    else return 0.0; 
  }
  else return vectorWedge2(pk, staticObsRep[alphaToAddOrUpdate(0)-num_of_agents-1].col(alphaToAddOrUpdate(1)), bp);  
}

double eu::getTetherLength(eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                         Eigen::Vector2d pb_self, Eigen::Vector2d& pkplus1, 
                         std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                         std::vector<Eigen::Vector2d>& staticObsLongestDist,
                         int num_of_agents)
{
  double length = 0.0;
  for (int i=0; i<entangle_state.bendPointsIdx.size(); i++)
  {
    Eigen::Vector2d bp;
    double compensation = 0;
    getBendPt2dwIdx(bp, entangle_state, pb, staticObsRep, staticObsLongestDist, num_of_agents, 
                    entangle_state.bendPointsIdx[i], compensation);
    length += (bp-pb_self).norm() + 2*compensation;
    pb_self = bp;
  }
  length += (pkplus1-pb_self).norm();

  return length;
}

double eu::tetherLength_orig(std::vector<Eigen::Vector2d>& contPts)
{
    double length = 0.0;
    for (int i=0; i<contPts.size()-1; i++){
        length += (contPts[i+1]-contPts[i]).norm();
    }
    return length;
}