/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkRectangleOptimizerv4_hxx
#define itkRectangleOptimizerv4_hxx

#include "itkRectangleOptimizerv4.h"
#include <cstdlib>
#include <time.h>
// #include </Users/cuicui/ITK/Modules/ThirdParty/VNL/src/vxl/core/vnl/algo/vnl_rectangle.h>

using namespace std;

float minf=100000;
float minarea=100000;
float minftemp=100000;
float minareatemp=100000;
int num=0;
class matrix m[1000];
float w0[3],cor0[3];
float maxgetrho;
int maxtemp=0;

float getfunc(float center[3])
{
    float f=0;
    return f;
}

void matrix::calculate(float w[3],float cor[3])
{
    int i;
    float GVn;
    for(i=0;i<3;i++)
    {
        width[i]=w[i];
        corner[i]=cor[i];
        center[i]=cor[i] + w[i]/2;
    }
    area=1;
    for(i=0;i<3;i++)
    {
        area*=width[i];
    }
    
    f=getfunc(center);
    
    if(num==0)
    {
        minf=f;
        minarea=area;
    }
    if(f<minf) minf=f;
    if(area<minarea) minarea=area;
    
    GVn = pow((minarea * log(1/minarea)), 2.0/3);
    getrho = pow(area,2.0/3)/(f-minf+3*GVn);
}

void matrix::update()
{
    float GVn;
    
    GVn = pow((minarea * log(1/minarea)), 2.0/3);
    getrho = pow(area,2.0/3)/(f-minf+3*GVn);
}

float * matrix::getwidth()
{
    float *temp=new float[3];
    int i;
    for(i=0;i<3;i++)
    {
        temp[i]=width[i];
    }
    return temp;
}
float * matrix::getcorner()
{
    float *temp=new float[3];
    int i;
    for(i=0;i<3;i++)
    {
        temp[i]=corner[i];
    }
    return temp;
}

float matrix::getgetrho()
{
    return getrho;
}

void matrix::savefile()
{
    FILE *f;
    int i;
    f=fopen("123.txt","a+");
    cout<<"center ";
    for(i=0;i<3;i++)
    {
        cout<<center[i]<<" ";
        fprintf(f,"%f ",center[i]);
    }
    cout<<endl;
    
    fprintf(f,"\n\n");
    fclose(f);
}


namespace itk
{

template<typename TInternalComputationValueType>
RectangleOptimizerv4<TInternalComputationValueType>
::RectangleOptimizerv4() :
  m_CurrentValue(0),
  m_NumberOfSteps(0),
  m_Stop(false),
  m_StepLength(0.5),
  m_CurrentIndex(0),
  m_MaximumMetricValue(0.0),
  m_MinimumMetricValue(0.0),
  m_StopConditionDescription("")
{
  this->m_NumberOfIterations = 0;
}

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::StartOptimization(bool /* doOnlyInitialization */)
{
  this->StartWalking();
}

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::StartWalking(void)
{
  itkDebugMacro("StartWalking");
  this->InvokeEvent( StartEvent() );
  m_StopConditionDescription.str("");
  m_StopConditionDescription << this->GetNameOfClass() << ": Running";
    
  ParametersType initialPos = this->m_Metric->GetParameters();
  m_MinimumMetricValuePosition = initialPos;
  m_MaximumMetricValuePosition = initialPos;

  this->SetInitialPosition( initialPos ); // store the initial position

  MeasureType initialValue = this->m_Metric->GetValue();
  m_MaximumMetricValue = initialValue;
  m_MinimumMetricValue = initialValue;
  
  this->m_CurrentIteration = 0;
  this->m_NumberOfIterations = 1;

  const unsigned int spaceDimension = this->m_Metric->GetParameters().GetSize();

  for ( unsigned int i = 0; i < spaceDimension; i++ )
    {
    this->m_NumberOfIterations *= ( 2 * m_NumberOfSteps[i] + 1 );
    }

  m_CurrentIndex.SetSize(spaceDimension);
  m_CurrentIndex.Fill(0);

  const ScalesType & scales = this->GetScales();
  // Make sure the scales have been set properly
  if ( scales.size() != spaceDimension )
    {
    itkExceptionMacro(<< "The size of Scales is "
                      << scales.size()
                      << ", but the NumberOfParameters is "
                      << spaceDimension
                      << ".");
    }

  // Setup first grid position.
    ParametersType position(spaceDimension);
 
    
    //cc add corner position
    //initialize width length and corner position
    //ParametersType corPos = this->m_Metric->GetParameters();
    
    for(int i=0;i<3;i++)
    {
        w0[i]=1;
        cor0[i]=0;
        position[i]=(cor0[i]+0.5*w0[i]-0.5)*scales[i]*m_NumberOfSteps[i] * m_StepLength;
    }
    m[0].calculate(w0,cor0);
    
    maxgetrho=m[0].getgetrho();
    maxtemp=0;
    
    this->m_Metric->SetParameters(position);
    m[num].f =this->m_Metric->GetValue();
    num++;
  itkDebugMacro("Calling ResumeWalking");

  this->ResumeWalking();
}
    
// the first rectangle
template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::ResumeWalking(void)
{
  itkDebugMacro("ResumeWalk");
  m_Stop = false;
    
  while ( !m_Stop )
    {
    ParametersType currentPosition = this->GetCurrentPosition();

    if ( m_Stop )
      {
      StopWalking();
      break;
      }

    m_CurrentValue = this->m_Metric->GetValue();

    if ( m_CurrentValue > m_MaximumMetricValue )
      {
      m_MaximumMetricValue = m_CurrentValue;
      m_MaximumMetricValuePosition = currentPosition;
      }
    if ( m_CurrentValue < m_MinimumMetricValue )
      {
      m_MinimumMetricValue = m_CurrentValue;
      m_MinimumMetricValuePosition = currentPosition;
      }
   //    const unsigned int spaceDimension = this->m_Metric->GetParameters().GetSize();

   //     ParametersType Width(spaceDimension);
   //     const ScalesType & scales = this->GetScales();


    if ( m_Stop )
      {
      this->StopWalking();
      break;
      }

//    m_StopConditionDescription.str("");
//    m_StopConditionDescription << this->GetNameOfClass() << ": Running. ";
//    m_StopConditionDescription << "@ index " << this->GetCurrentIndex() << " value is " << m_CurrentValue;

    this->InvokeEvent( IterationEvent() );
  // this->AdvanceOneStep();
         this->Advance();
    this->m_CurrentIteration++;
//       if(this->m_CurrentIteration>530)
//        {
//            m_Stop = true;
//        }
    }
}

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::StopWalking(void)
{
  itkDebugMacro("StopWalking");

  m_Stop = true;
  this->InvokeEvent( EndEvent() );
}
    

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::Advance(void)
{
    
    int j,k;
    
    float *maxwidth,*maxcorner;
    float tempwidth[3][3],tempcorner[3][3];
        itkDebugMacro("Advance");
        
        const unsigned int spaceDimension = this->m_Metric->GetParameters().GetSize();
        
        ParametersType newPosition(spaceDimension);
        //this->Increment(newPosition);
           float newposition[3];
    ParametersType position1(spaceDimension);
     ParametersType position2(spaceDimension);
     ParametersType position3(spaceDimension);
    maxwidth=m[maxtemp].getwidth();
    maxcorner=m[maxtemp].getcorner();
    
    //slipt one rectangle into three rectangles
    float tempmaxwidth=0;
    int tempmax;
    for(j=0;j<3;j++)
    {
        if(maxwidth[j]>tempmaxwidth)//slipt by largest width
        {
            tempmaxwidth=maxwidth[j];
            tempmax=j;
        }
    }
    
    for(j=0;j<3;j++)
    {
        tempwidth[0][j]=maxwidth[j];
        tempwidth[1][j]=maxwidth[j];
        tempwidth[2][j]=maxwidth[j];
    }
    tempwidth[0][tempmax]=maxwidth[tempmax]/3;
    tempwidth[1][tempmax]=maxwidth[tempmax]/3;
    tempwidth[2][tempmax]=maxwidth[tempmax]/3;
    
    for(j=0;j<3;j++)
    {
        tempcorner[0][j]=maxcorner[j];
        tempcorner[1][j]=maxcorner[j];
        tempcorner[2][j]=maxcorner[j];
    }
    tempcorner[0][tempmax]=maxcorner[tempmax];
    tempcorner[1][tempmax]=maxcorner[tempmax]+maxwidth[tempmax]/3;
    tempcorner[2][tempmax]=maxcorner[tempmax]+maxwidth[tempmax]/3*2;
    
    
    minftemp=minf;
    minareatemp=minarea;
    
    
    m[num].calculate(tempwidth[0],tempcorner[0]);
    for(k=0;k<3;k++)
    {
        position1[k]=m[num].center[k];
    }
    this->m_Metric->SetParameters(position1);
    m[num].f =this->m_Metric->GetValue();
  
    num++;
   
    m[num].calculate(tempwidth[1],tempcorner[1]);
    for(k=0;k<3;k++)
    {
        position2[k]=m[num].center[k];
    }
    this->m_Metric->SetParameters(position2);
    m[num].f =this->m_Metric->GetValue();
 
    num++;
    
    m[num].calculate(tempwidth[2],tempcorner[2]);
    for(k=0;k<3;k++)
    {
        position3[k]=m[num].center[k];
    }
    this->m_Metric->SetParameters(position3);
    m[num].f =this->m_Metric->GetValue();
 
    num++;
    
    
    for(j=maxtemp;j<num;j++)
    {
        m[j]=m[j+1];
    }
    num--;
    
    //update rectangles
    if(minftemp != minf || minareatemp != minarea)
    {
        for(j=0;j<num;j++)
        {
            m[j].update();
        }
    }
    
    //getrho
    maxgetrho=0;
    if(num==1)
    {
        maxgetrho=m[0].getgetrho();
        maxtemp=0;
    }
    else
    {
        for(j=0;j<num;j++)
        {
//              cout<<"num: "<<j+1<<" getrho:"<<m[j].getgetrho()<<" getf:"<<m[j].f<<" ";
//             for(k=0;k<3;k++)
//             {
//             cout<<m[j].center[k]<<" ";
//             }
//             cout<<endl;
            
            if(m[j].getgetrho()>maxgetrho)
            {
                maxgetrho=m[j].getgetrho();
                maxtemp=j;//the one who has the largest rho
                for(k=0;k<3;k++)
                {
                    newposition[k]=m[j].center[k];
                }
            }
        }
    }
        const ScalesType & scales = this->GetScales();
    
        for ( unsigned int i = 0; i < spaceDimension; i++ )
        {
            newPosition[i]= (newposition[i]-0.5)*scales[i]*m_NumberOfSteps[i] * m_StepLength;
        //    newPosition[i]=newposition[i];
            
        }
    
    
        //  newPosition = newposition;
        itkDebugMacro(<< "new position = " << newPosition);
        
        this->m_Metric->SetParameters(newPosition);
        
}

template<typename TInternalComputationValueType>
const std::string
RectangleOptimizerv4<TInternalComputationValueType>
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::SetInitialPosition(const ParametersType & param)
{
  m_InitialPosition = param;
  this->Modified();
}

template<typename TInternalComputationValueType>
void
RectangleOptimizerv4<TInternalComputationValueType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "InitialPosition = " << m_InitialPosition << std::endl;
  os << indent << "CurrentValue = " << m_CurrentValue << std::endl;
  os << indent << "NumberOfSteps = " << m_NumberOfSteps << std::endl;
  os << indent << "Stop = " << m_Stop << std::endl;
  os << indent << "StepLength = " << m_StepLength << std::endl;
  os << indent << "CurrentIndex = " << m_CurrentIndex << std::endl;
  os << indent << "MaximumMetricValue = " << m_MaximumMetricValue << std::endl;
  os << indent << "MinimumMetricValue = " << m_MinimumMetricValue << std::endl;
  os << indent << "MinimumMetricValuePosition = " << m_MinimumMetricValuePosition << std::endl;
  os << indent << "MaximumMetricValuePosition = " << m_MaximumMetricValuePosition << std::endl;
}
 
} // end namespace itk

#endif
