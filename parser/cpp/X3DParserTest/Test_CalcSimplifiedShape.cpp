#include "CX3DParser.h"
#include "CSimplifiedShape.h"
#include "Test.h"
#include <vector>
#include <string>

// ============================================
// ============================================

void Test_calcSimplifiedShape(CX3DParser& parser, char *vrmlFile)
{
	// ---------------------------------------------
	// ---------------------------------------------
	if (!parser.parse(vrmlFile))
	{
		fprintf(stderr, "%s parse failed\n", vrmlFile);
		return;
	}

	// ---------------------------------------------
	//
	//
	//  Group
	//    Transform
	//      Shape
	//        IndexedFaceSet
	// ---------------------------------------------

	MFNode *groupNodes = parser.searchNodesFromDirectChildrenOfRoot("Group");
	if (!groupNodes)
	{
		fprintf(stderr, "could not found root group node in %s\n", vrmlFile);
		return;
	}

	CX3DGroupNode *pGroup = (CX3DGroupNode *)(groupNodes->getNode(0));
	if (!pGroup)
	{
		fprintf(stderr, "root group is NULL\n");
		return;
	}

	MFNode *transformNodes = pGroup->searchNodesFromDirectChildren("Transform");
	if (!transformNodes)
	{
		fprintf(stderr, "could not found Transform node\n");
		return;
	}

	CX3DTransformNode *pTrans = (CX3DTransformNode *)(transformNodes->getNode(0));
	if (!pTrans)
	{
		fprintf(stderr, "Transform node is NULL\n");
		return;
	}

	MFNode *shapeNodes = pTrans->searchNodesFromDirectChildren("Shape");
	if (!shapeNodes)
	{
		fprintf(stderr, "could not found Shape node\n");
		return;
	}

	// -------------------------------------------
	// -------------------------------------------

	CSimplifiedShape *ss = CSimplifiedShapeFactory::calc(shapeNodes);
	if (!ss)
	{
		fprintf(stderr, "error in calculating simplified shape\n");
		return;
	}

	switch (ss->getType())
	{
	case CSimplifiedShape::SPHERE:
		{
			CSimplifiedSphere *sp = (CSimplifiedSphere *)ss;

			printf("** sphere **\n");
			printf("\tradius  : %f\n", sp->radius());
			printf("\tcenter  : (%f, %f, %f)\n", sp->x(), sp->y(), sp->z());
			printf("\t�c�x : %f\n", sp->hizumi());
		}
		break;

	case CSimplifiedShape::CYLINDER:
		{
			CSimplifiedCylinder *cy = (CSimplifiedCylinder *)ss;

			printf("** cylinder **\n");
			printf("\tradius  : %f\n", cy->radius());
			printf("\theight  : %f\n", cy->height());
			printf("\tcenter  : (%f, %f, %f)\n", cy->x(), cy->y(), cy->z());
			printf("\t�c�x : %f\n", cy->hizumi());
		}
		break;
		
	case CSimplifiedShape::BOX:
		{
			CSimplifiedBox *bx = (CSimplifiedBox *)ss;

			printf("** box **\n");
			printf("\tp1  : (%f, %f, %f)\n", bx->x1(), bx->y1(), bx->z1());
			printf("\tp2  : (%f, %f, %f)\n", bx->x2(), bx->y2(), bx->z2());
			printf("\t�c�x : %f\n", bx->hizumi());
		}
		break;
	}

	delete ss;


	printf("�e�P���`��̘c�x\n");
	CSimplifiedSphere *sp = CSimplifiedShapeFactory::calcSphere(shapeNodes);
	if (sp)
	{
		printf("\t��         : %f\n", sp->hizumi());

		delete sp;
	}

	CSimplifiedCylinder *cy = CSimplifiedShapeFactory::calcCylinder(shapeNodes);
	if (cy)
	{
		printf("\t�V�����_�[ : %f\n", cy->hizumi());

		delete cy;
	}

	CSimplifiedBox *bx = CSimplifiedShapeFactory::calcBox(shapeNodes);
	if (bx)
	{
		printf("\t��         : %f\n", bx->hizumi());

		delete bx;
	}

	delete shapeNodes;
	delete transformNodes;
	delete groupNodes;
}

