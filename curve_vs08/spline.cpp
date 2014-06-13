#include "spline.h"
#include <vector>

const int Division = 20;

Spline::Spline(GraphWidget *graphWidget)
:graph(graphWidget)
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);

	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);

	startPoint = NULL;
	endPoint = NULL;
	totalPoints = 0;
	style = BezierBerstein;
	showControlPoints = false;
	showPseudoPoints = true;
}

void Spline::SetVisible(QList<Node *> nodes, bool visible)
{
	foreach(Node * i, nodes)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(QList<Edge *> edges, bool visible)
{
	foreach(Edge * i, edges)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(bool visible)
{
	if (this->startPoint && this->endPoint)
	{
		this->startPoint->setVisible(visible);
		this->startPoint->setEnabled(visible);
		this->endPoint->setVisible(visible);
		this->endPoint->setEnabled(visible);
		if (!visible)
		{
			this->graph->scene()->removeItem(startPoint);
			this->graph->scene()->removeItem(endPoint);
		}
		else
		{
			this->graph->scene()->addItem(startPoint);
			this->graph->scene()->addItem(endPoint);
		}
	}
}

void Spline::Clear()
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);

	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);

	if (startPoint)
	{
		this->graph->scene()->removeItem(startPoint);
		delete startPoint;
		startPoint = NULL;
	}
	if (endPoint)
	{
		this->graph->scene()->removeItem(endPoint);
		delete endPoint;
		endPoint = NULL;
	}

	totalPoints = 0;
}

void Spline::PaintPseudeLines()
{
	if (startPoint && endPoint)
	{
		Interpolate();
		deleteEdges(pseudoEdges);
		calcPseudoEdges();
		Paint();
	}
}

void Spline::deleteNodes(QList<Node *> &nodes)
{
	foreach(Node * i, nodes)
	{
		delete i;
	}
	nodes.clear();
}

void Spline::deleteEdges(QList<Edge *> &edges)
{
	foreach(Edge * i, edges)
	{
		delete i;
	}
	edges.clear();
}

void Spline::Paint()
{
	deleteEdges(curveEdges);
	deleteEdges(pseudoEdges);
	deleteEdges(controlEdges);

	// Curve Points & Edges
	foreach(Node * i, curvePoints)
	{
		this->graph->scene()->addItem(i);
	}
	calcCurveEdges();

	// Control Points & Edges
	if (showControlPoints)
	{
		foreach(Node * i, controlPoints)
		{
			this->graph->scene()->addItem(i);
		}
		calcControlEdges();
	}

	if (showPseudoPoints)
	{
		calcPseudoEdges();
		this->SetVisible(true);
	}
	else
	{
		this->SetVisible(false);
	}

	this->graph->repaint(-1000, -1000, 2000, 2000);
}

void Spline::AddPoint(Node *newNode)
{
	interpPoints.append(newNode);
	graph->scene()->addItem(newNode);
	totalPoints++;

	CalcStartPoint();
	CalcEndPoint();
}

void Spline::CalcStartPoint()
{
	if (totalPoints == 2)
	{
		Node* p1 = interpPoints[0];
		Node* p2 = interpPoints[1];

		QPointF line12 = p1->pos() - p2->pos();
		double length12 = sqrt(line12.x()*line12.x() + line12.y()*line12.y());

		startPoint = new Node(this->graph, PseudoPoints);
		startPoint->setPos(p1->pos() + line12 / length12 * 50);
		graph->scene()->addItem(startPoint);
	}
}

void Spline::CalcEndPoint()
{
	if (totalPoints < 2)
		return;
	if (!endPoint)
	{
		endPoint = new Node(this->graph, PseudoPoints);
		graph->scene()->addItem(endPoint);
	}

	Node* pN = interpPoints[totalPoints - 1];
	Node* pNm1 = interpPoints[totalPoints - 2];

	QPointF lineNNm1 = pN->pos() - pNm1->pos();
	double lengthNNm1 = sqrt(lineNNm1.x()*lineNNm1.x() + lineNNm1.y()*lineNNm1.y());

	endPoint->setPos(pN->pos() + lineNNm1 / lengthNNm1 * 50);
}

void Spline::calcCurveEdges()
{
	deleteEdges(curveEdges);
	for (int i = 0; i < curvePoints.size() - 1; i++)
	{
		Edge * tmp = new Edge(curvePoints[i]->pos(), curvePoints[i+1]->pos());
		curveEdges.append(tmp);
		graph->scene()->addItem(tmp);
	}
}

void Spline::calcControlEdges()
{
	deleteEdges(controlEdges);

	switch(style)
	{
	case BezierBerstein:
	case BezierDeCasteljau:
	case BezierMatrix:
		if (this->controlPoints.size() == 2*(totalPoints-1))
		{
			for (int i = 0; i < this->totalPoints - 1; i ++)
			{
				Node * n1 = interpPoints[i];
				Node * n2 = interpPoints[i+1];
				Node * c1 = controlPoints[2*i];
				Node * c2 = controlPoints[2*i+1];
				Edge * e1 = new Edge(n1->pos(), c1->pos(), 1);
				Edge * e2 = new Edge(c1->pos(), c2->pos(), 1);
				Edge * e3 = new Edge(c2->pos(), n2->pos(), 1);
				controlEdges.append(e1);
				controlEdges.append(e2);
				controlEdges.append(e3);
				graph->scene()->addItem(e1);
				graph->scene()->addItem(e2);
				graph->scene()->addItem(e3);
			}
		}
		break;
	case BSpline:	
		for (int i = 0; i < this->controlPoints.size() - 1; i++)
		{
			Node * c1 = controlPoints[i];
			Node * c2 = controlPoints[i+1];

			Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
			controlEdges.append(e1);
			graph->scene()->addItem(e1);
		}
		break;
	case HermiteClamped:
	case HermiteNatural:
		if (this->controlPoints.size() == totalPoints)
		{
			for (int i = 0; i < this->totalPoints; i++)
			{
				Node * c1 = interpPoints[i];
				Node * c2 = controlPoints[i];

				c2->setPosition(c1->getPosition() + c2->getPosition());

				Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
				controlEdges.append(e1);
				graph->scene()->addItem(e1);
			}
		}
		break;
	}
} 

void Spline::calcPseudoEdges()
{
	if (totalPoints >= 2)
	{
		deleteEdges(pseudoEdges);
		// start
		Edge * start = new Edge(startPoint->pos(), ((Node *)interpPoints.first())->pos(), 2);
		pseudoEdges.append(start);
		graph->scene()->addItem(start);
		// end
		Edge * end = new Edge(endPoint->pos(), ((Node *)interpPoints.last())->pos(), 2);
		pseudoEdges.append(end);
		graph->scene()->addItem(end);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Start your code here.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Spline::Interpolate()
{
	//Clear the old curve points
	deleteNodes(curvePoints);
	//Calculate control points.
	CalcControlPoints();
	//Depending on the selected style, interpolate the curve.
	switch(style){
	case BezierBerstein:		InterpBerstein();       break;
	case BezierDeCasteljau:		InterpCasteljau();      break;
	case BezierMatrix:          InterpMatrix();         break;
	case BSpline:               InterpBSpline();        break;
	case HermiteClamped:        InterpHermiteClamped(); break;
	case HermiteNatural:        InterpHermiteNatural(); break;
	}  
}

float Spline::N(int n, int j, float t)
{  
	vector<float> vLamda (1000); 
	//store all the points into lamda vector
	vLamda[0] = -3; 
	vLamda[1] = -2; 
	vLamda[2] = -1; 

	//think about the limit here 
	for(int i=0; i< totalPoints; i++)
	{
		vLamda[i+3] = i; 
	}
	vLamda[totalPoints+3] = totalPoints;
	vLamda[totalPoints+4] = totalPoints+1;
	vLamda[totalPoints+5] = totalPoints+2; 

	//there should be a comparison about j here, am I correct? 
	if((n==0)&&(t>=vLamda[j])&&(t<vLamda[j+1]))
	{
		return 1.0; 
	}
	else if((n==0)&&((t<vLamda[j])||(t>=vLamda[j+1])))
	{
		return 0.0; 
	}
	else if((t<vLamda[j])||(t>vLamda[j+n+1]))
	{
		return 0.0; 
	}else {
		float firstPart = ((float)(t-vLamda[j])/(float)(vLamda[j+n]-vLamda[j]))*N(n-1, j, t);
		float secondPart = ((float)(vLamda[j+n+1]-t)/(float)(vLamda[j+n+1]-vLamda[j+1]))*N(n-1, j+1, t);
		float final = firstPart+ secondPart; 
		return final; 
	}
}  

float Spline::Nder(int l, int n, int j, float t)
{
	vector<float> vLamda (1000); 
	//store all the points into lamda vector
	vLamda[0] = -3; 
	vLamda[1] = -2; 
	vLamda[2] = -1; 

	//think about the limit here 
	for(int i=0; i< totalPoints; i++)
	{
		vLamda[i+3] = i; 
	}
	vLamda[totalPoints+3] = totalPoints;
	vLamda[totalPoints+4] = totalPoints+1;
	vLamda[totalPoints+5] = totalPoints+2; 

	if(l==0)
	{
		return N(n, j, t); 
	}
	else{
		float firstPart = ((float)1.0/(vLamda[j+n]-vLamda[j]))*Nder(l-1,n-1,j,t);
		float secondPart = ((float)1.0/(vLamda[j+n+1]-vLamda[j+1]))*Nder(l-1, n-1, j+1, t);
		float final = (float) (n*(firstPart-secondPart)); 
		return final;  
	}
}

//////////////////////////////////////////////////////////////////////////
// Calculate the control points
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// controlPoints- type: QList<Node *>
//				  discription: stores all the control points for curve interpolation and display.
//                             For Bezier curve, between very pair of consecutive interpolation points,
//                             there should be two control points. These four points determins the curve interpolation.
//                             For B-Spline, there should be totalPoints + 2 control points calculated from Ac = p.
//
// Hint: To implement Hermite and B-Spline, you need to write functions to create the A matrix as in the handouts.
//       Then you solve a linear system Ac = p, where p is the interpolation points vector and c are the control points.
//       We have provided you with a data structure to store and solve the linear system.
//       Below is an example code, read the understand it.
//
//	matrix<float> A(3,3);
//  matrix<float> c(3,1);
//  matrix<float> p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
//  c = A.Solve(p);
//  
//  The result in c is c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.
// Hint2: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 startPointPos = startPoint->getPosition();
// 
// Hint3: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example: 
// Node* ctrl = new Node(graph, ControlPoints);
// ctrl->setPosition(/*some vec3*/);
// controlPoints.append(ctrl);
void Spline::CalcControlPoints()
{
	// Clear previous controlPoints
	deleteNodes(controlPoints);
	// Based on style use appropriate ways to compute control points
	if(totalPoints<=1)
	{
		return; 
	}

	switch(style)
	{
	case BezierBerstein:
	case BezierDeCasteljau:
	case BezierMatrix:
	{
			// ADD YOUR CODE HERE
			vec3 startPointPos = startPoint->getPosition(); 
			vec3 endPointPos = endPoint->getPosition(); 
			//set the first controlPoint
			Node* ctrl = new Node(graph, ControlPoints); 
			vec3 control0 = interpPoints[0]->getPosition() + (interpPoints[0]->getPosition() - startPointPos); 
			ctrl->setPosition(control0); 
			controlPoints.append(ctrl); 
			for (int i = 1; i <= this->totalPoints - 2; i++)
			{
				//control[2i-1]
				Node* ctrl1 = new Node(graph, ControlPoints); 
				vec3 contro2i1 = interpPoints[i]->getPosition()- ((interpPoints[i+1]->getPosition() - interpPoints[i-1]->getPosition())/6);
				ctrl1->setPosition(contro2i1);
				controlPoints.append(ctrl1); 
				//control[2i]
				Node* ctrl2 = new Node(graph, ControlPoints); 
				vec3 contro2i =  interpPoints[i]->getPosition() + ((interpPoints[i+1]->getPosition() -interpPoints[i-1]->getPosition())/6); 
				ctrl2->setPosition(contro2i); 
				controlPoints.append(ctrl2);  
			} 
			Node* ctrl3 = new Node(graph, ControlPoints);
			vec3 controlast = interpPoints[totalPoints-1]->getPosition() + (interpPoints[totalPoints-1]->getPosition() - endPointPos);
			ctrl3->setPosition(controlast); 
			controlPoints.append(ctrl3); 
			break; 
	}
		break;
	case BSpline:
		// make sure there are at least 2 totalPoints.
		// ADD YOUR CODE HERE  
		{
		matrix<float> A(totalPoints+2, totalPoints+2); 
		matrix<float> c(totalPoints+2, 2); 
		matrix<float> p(totalPoints+2, 2); 
		A(0,0) = (float)(Nder(2, 3, 0, 0));
		A(0,1) = (float)(Nder(2, 3, 1, 0));
		A(0,2) = (float)(Nder(2, 3, 2, 0));
		for(int i =3; i<totalPoints+2; i++)
		{      
				A(0,i)= 0; 
		}
		 for(int i=1; i<totalPoints+1; i++)
		 {
			   for(int j=0; j<totalPoints+2; j++)
			   {
				   A(i,j) = N(3, j, i-1);
			   }
		 }   
		 //store the last line
		 for(int i=0; i<totalPoints-1; i++)
		 {
				 A(totalPoints+1, i) = 0; 
		 }
		 A(totalPoints+1, totalPoints-1) =  (float)(Nder(2, 3, totalPoints-1, totalPoints-1)); 
		 A(totalPoints+1, totalPoints) =  (float)(Nder(2, 3, totalPoints, totalPoints-1)); 
		 A(totalPoints+1, totalPoints+1) = (float)(Nder(2, 3, totalPoints+1, totalPoints-1)); 
		 p(0,0) = 0.0; 
		 p(0,1) = 0.0; 
		 p(totalPoints+1, 0)=0.0;
		 p(totalPoints+1, 1)=0.0;
		 for(int i=1; i<totalPoints+1; i++)
		 {
			   p(i,0) = interpPoints[i-1]->x();
			   p(i,1) = interpPoints[i-1]->y(); 
		 }
		 c = A.Solve(p);
		 for(int i=1; i<totalPoints+1; i++)
		 {
			    vec3 cValue(c(i,0), c(i,1), 0); 
			    Node* ctrl3 = new Node(graph, ControlPoints);
			    ctrl3->setPosition(cValue); 
			    controlPoints.append(ctrl3); 
		 } 
		break;
	}
	// In the case for Hermite, you want to implement both "clamped" and "natural" versions. 
	case HermiteClamped:
		// make sure there are at least 2 totalPoints.
		// ADD YOUR CODE HERE
	{
		matrix<float> A(totalPoints, totalPoints); 
		matrix<float> c(totalPoints, 2); 
		matrix<float> p(totalPoints, 2); 
		float startPointX = startPoint->x(); 
		float startPointY = startPoint->y(); 
		float endPointX = endPoint->x(); 
		float endPointY = endPoint->y(); 

		   //store the first line 
		   A(0,0) = 1; 
		   for(int i =1; i<totalPoints; i++)
		   {
			  A(0,i)=0;
		   }
		   //store the middle parts
		   for(int i=1; i<totalPoints-1; i++)
		   {
			   for(int j=0; j<totalPoints; j++)
			   {
				   if(i==(j+1))
				   {
                      A(i,j) = 1; 
				   }else if(i==j)
				   { 
                      A(i,j) = 4; 
				   }else if(i==(j-1))
				   {
					  A(i,j) = 1; 
				   }else{
					  A(i,j) = 0; 
				   }
			   }
		   } 
		   //store the last line
		   for(int i=0; i<totalPoints-1; i++)
		   {
			   A(totalPoints-1, i) = 0;
		   }
		   A(totalPoints-1, totalPoints-1) = 1; 

		   p(0,0) = startPointX - interpPoints[0]->x(); 
		   p(0,1) = startPointY - interpPoints[0]->y(); 
		   p(totalPoints-1, 0)=endPointX - interpPoints[totalPoints-1]->x();
		   p(totalPoints-1, 1)=endPointY - interpPoints[totalPoints-1]->y();
		   for(int i=1; i<totalPoints-1; i++)
		   {
			   p(i,0) = 3*(interpPoints[i+1]->x() - interpPoints[i-1]->x());
			   p(i,1) = 3*(interpPoints[i+1]->y() - interpPoints[i-1]->y()); 
		   }
		   //solve the problem 
		   c = A.Solve(p);
		   for(int i=0; i<totalPoints; i++)
		   {
			    vec3 cValue(c(i,0), c(i,1), 0); 
			    Node* ctrl3 = new Node(graph, ControlPoints);
			    ctrl3->setPosition(cValue); 
			    controlPoints.append(ctrl3); 
		   } 
		break;
	}
	
	case HermiteNatural:
		{
		// make sure there are at least 2 totalPoints.
		// ADD YOUR CODE HERE
		matrix<float> A(totalPoints, totalPoints); 
		matrix<float> c(totalPoints, 2); 
		matrix<float> p(totalPoints, 2); 
		float startPointX = startPoint->x(); 
		float startPointY = startPoint->y(); 
		float endPointX = endPoint->x(); 
		float endPointY = endPoint->y();
		//store a
		A(0,0) = 2; A(0,1) = 1; 
	    for(int i =2; i<totalPoints; i++)
		{
			  A(0,i)=0;
	    }
		 for(int i=1; i<totalPoints-1; i++)
		 {
			   for(int j=0; j<totalPoints; j++)
			   {
				   if(i==(j+1))
				   {
                      A(i,j) = 1; 
				   }else if(i==j)
				   { 
                      A(i,j) = 4; 
				   }else if(i==(j-1))
				   {
					  A(i,j) = 1; 
				   }else{
					  A(i,j) = 0; 
				   }
			   }
		 }   
		 //store the last line
		 for(int i=0; i<totalPoints-2; i++)
		 {
			   A(totalPoints-1, i) = 0;
		 }
		 A(totalPoints-1, totalPoints-2) = 1; 
		 A(totalPoints-1, totalPoints-1) = 2; 

		 p(0,0) = 3*(interpPoints[1]->x() - interpPoints[0]->x()); 
		 p(0,1) = 3*(interpPoints[1]->y() - interpPoints[0]->y()); 
		 p(totalPoints-1, 0)=3*(interpPoints[totalPoints-1]->x() - interpPoints[totalPoints-2]->x());
		 p(totalPoints-1, 1)=3*(interpPoints[totalPoints-1]->y() - interpPoints[totalPoints-2]->y());
		 for(int i=1; i<totalPoints-1; i++)
		 {
			   p(i,0) = 3*(interpPoints[i+1]->x() - interpPoints[i-1]->x());
			   p(i,1) = 3*(interpPoints[i+1]->y() - interpPoints[i-1]->y()); 
		 }
		 c = A.Solve(p);
		 for(int i=0; i<totalPoints; i++)
		 {
			vec3 cValue(c(i,0), c(i,1), 0); 
			Node* ctrl3 = new Node(graph, ControlPoints);
			ctrl3->setPosition(cValue); 
			controlPoints.append(ctrl3); 
		 }
		break; 
		}
    }
}

//////////////////////////////////////////////////////////////////////////
// Cubic Berstein Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
// 
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
void Spline::InterpBerstein()
{
	// ADD YOUR CODE HERE
	if(totalPoints<=1)
	{
		return; 
	}
	for(int i= 0; i<=totalPoints-2; i++)
	{
		for(float u = 0.0; u<1.00; u=u+0.05)
		{
			float B0 = (1-u)*(1-u)*(1-u); 
			float B1 = 3*u*(1-u)*(1-u); 
			float B2 = 3*u*u*(1-u); 
			float B3 = u*u*u; 
			vec3 Qu = B0*interpPoints[i]->getPosition()+B1*controlPoints[2*i]->getPosition()+B2*controlPoints[2*i+1]->getPosition()+B3*interpPoints[i+1]->getPosition();
			Node* newNode = new Node(graph, CurvePoints); 
			newNode->setPosition(Qu);
			curvePoints.append(newNode);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Cubic de Casteljau Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
// 
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
void Spline::InterpCasteljau()
{
	// ADD YOU CODE HERE
	for(int i = 0; i<=totalPoints-2; i++)
	{
		for(float u=0.0; u<1.00; u=u+0.05)
		{
			vec3 P00 = interpPoints[i]->getPosition();
			vec3 P01 = controlPoints[2*i]->getPosition(); 
			vec3 P02 = controlPoints[2*i+1]->getPosition();
			vec3 P03 = interpPoints[i+1]->getPosition(); 
			
			vec3 P10 = (1-u)*P00+u*P01;
			vec3 P11 = (1-u)*P01+u*P02; 
			vec3 P12 = (1-u)*P02+u*P03; 

			vec3 P20 = (1-u)*P10+u*P11; 
			vec3 P21 = (1-u)*P11+u*P12; 
			
			vec3 P30 = (1-u)*P20+u*P21; 

			vec3 Qu = P30; 

			Node* newNode = new Node(graph, CurvePoints); 
			newNode->setPosition(Qu);
			curvePoints.append(newNode);
		}
	} 
} 

//////////////////////////////////////////////////////////////////////////
// Cubic de Casteljau Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
// 
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
// Hint3: If you have any problem in vec3 vec4 mat4 operations, see the SvzAlgebra.h and svzalgebra.cpp to find some specific method out.
// For Example:
// vec3 operator * (const mat4& a, const vec3& v); has been defined in SvzAlgebra.h
// That means you can do operations like this:
// vec3 A;
// mat4 M;
// vec3 B = M*A;
// (It automatically transfer the A and B vectors to homogeneous vectors)
void Spline::InterpMatrix()
{
	// ADD YOUR CODE HERE.
	if(totalPoints<=1)
	{
		return; 
	}
	for(int i= 0; i<=totalPoints-2; i++)
	{
		for(float u = 0.0; u<1.00; u=u+0.05)
		{
			vec4 firstRow(-1, 3, -3, 1);
			vec4 secondRow(3, -6, 3, 0);
			vec4 thirdRow(-3,3,0,0); 
			vec4 fourthRow(1, 0, 0, 0); 
			mat4 M(firstRow, secondRow, thirdRow, fourthRow); 
			vec4 U(u*u*u, u*u, u, 1.0);
			mat4 G(interpPoints[i]->getPosition(), controlPoints[2*i]->getPosition(), controlPoints[2*i+1]->getPosition(), interpPoints[i+1]->getPosition());
			
			vec3 Qu = (U*M)*G; 
			Node* newNode = new Node(graph, CurvePoints); 
			newNode->setPosition(Qu);
			curvePoints.append(newNode);
			 
		}
	}
    
}

//////////////////////////////////////////////////////////////////////////
// BSpline curve
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
void Spline::InterpBSpline()
{
	// ADD YOUR CODE HERE
	if(totalPoints<=1)
	{
		return; 
	}

	for(float t=0.0; t<=(totalPoints-1); t=t+0.05)
	{
		int integerT = (int) (t); 
		int j = integerT+3; 
		float N0 = N(3, j-3, t); 
		float N1 = N(3, j-2, t); 
		float N2 = N(3, j-1, t); 
		float N3 = N(3, j, t);
		vec3 Ft;
		if(t<(totalPoints-4))
		{
		    Ft = N0*controlPoints[j-3]->getPosition()+N1*controlPoints[j-2]->getPosition()+N2*controlPoints[j-1]->getPosition()+N3*controlPoints[j]->getPosition();
		}
	
		Node* newNode = new Node(graph, CurvePoints); 
	    newNode->setPosition(Ft);
		curvePoints.append(newNode);
	}

}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Clamped' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
void Spline::InterpHermiteClamped()
{
	// ADD YOUR CODE HERE
	if(totalPoints<=1)
	{
		return; 
	}
	for(int i= 0; i<=totalPoints-2; i++)
	{
		for(float u = 0.0; u<1.00; u=u+0.05)
		{
		    float H0 = 2*u*u*u - 3*u*u + 1; 
			float H1 = -2*u*u*u + 3*u*u; 
			float H2 = u*u*u-2*u*u+u; 
			float H3 = u*u*u-u*u; 
			vec3 Qu = H0*interpPoints[i]->getPosition()+ H1*interpPoints[i+1]->getPosition() +H2*controlPoints[i]->getPosition()+H3*controlPoints[i+1]->getPosition();
			Node* newNode = new Node(graph, CurvePoints); 
			newNode->setPosition(Qu);
			curvePoints.append(newNode);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Natural' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables 
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
// Hint: It is exactly the same like the method above.
void Spline::InterpHermiteNatural()
{
	// ADD YOUR CODE HERE.
	if(totalPoints<=1)
	{
		return; 
	}
	for(int i= 0; i<=totalPoints-2; i++)
	{
		for(float u = 0.0; u<1.00; u=u+0.05)
		{
		    float H0 = 2*u*u*u - 3*u*u + 1; 
			float H1 = -2*u*u*u + 3*u*u; 
			float H2 = u*u*u-2*u*u+u; 
			float H3 = u*u*u-u*u; 
			vec3 Qu = H0*interpPoints[i]->getPosition()+ H1*interpPoints[i+1]->getPosition() +H2*controlPoints[i]->getPosition()+H3*controlPoints[i+1]->getPosition();
			Node* newNode = new Node(graph, CurvePoints); 
			newNode->setPosition(Qu);
			curvePoints.append(newNode);
		}
	}
}
