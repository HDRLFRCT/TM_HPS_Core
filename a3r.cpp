//-----------------------------------------------------------------------------------
//
//	Saving and loading particle coordinates in A3R format
//
//	Anton M. Krivtsov
//
//	29.10.2001
//
//-------------------------------------------------------------------------------------------

#include "a3r.h"
#include <fstream>

//-------------------------------------------------------------------------------------------

A3R_HEADER::A3R_HEADER()
: count(0)
, r(0)
, count_1(0)
{
	strcpy(file_type, "a3r");
	strcpy(version, "a");
	data_start = sizeof(A3R_HEADER);
}

//-------------------------------------------------------------------------------------------

BOOL Save_A3R(const char* file_name, Vect3D* start, LONG n, FLT r)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::out)) return FALSE;
#else
	if(!file.open(file_name, ios::out | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	header.count = n;
	header.r = r;

	file.sputn((char*) &header, sizeof(header));
	for (Vect3D* i = start; i != start + n; i++)
	{
		Vect3D::SFLOAT sr = *i;
		file.sputn((char*) &sr, sizeof(Vect3D::SFLOAT));
	}

	file.close();

	return TRUE;
}

//-------------------------------------------------------------------------------------------

Vect3D* Load_A3R(const char* file_name, LONG& n, FLT& r)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::in | ios::nocreate)) return NULL;
#else
//	if(!file.open(file_name, ios::in | ios::nocreate | ios::binary)) return NULL;
#endif

	A3R_HEADER header;

	file.sgetn((char*)&header, sizeof(A3R_HEADER));
	if(strcmp(header.file_type, "a3r")) return NULL;
	n = header.count;
	r = float(header.r);

	Vect3D::SFLOAT* buf = new Vect3D::SFLOAT[n];
	Vect3D::SFLOAT* j = buf;

	file.sgetn((char*)buf, sizeof(Vect3D::SFLOAT) * n);

	file.close();

	Vect3D* start = new Vect3D[n];
	Vect3D* stop = start + n;
	for (Vect3D* i = start; i != stop; i++) *i = *j++;
	
	delete [] buf;

	return start;
}



//-------------------------------------------------------------------------------------------

BOOL Save_A3R(const char* file_name, Atom2D* as, LONG n, FLT r)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::out)) return FALSE;
#else
	if(!file.open(file_name, ios::out | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	header.count = n;
	header.r = r;

	file.sputn((char*) &header, sizeof(header));
	for (Atom2D* i = as; i != as + n; i++)
	{
		Vect3D r(i->r.x, i->r.y, 0);
	//	r.RotateX(DEG * .001);
	//	r.RotateY(DEG * .001);
		Vect3D::SFLOAT sr = r;
		file.sputn((char*) &sr, sizeof(Vect3D::SFLOAT));
	}

	file.close();

	return TRUE;
}

BOOL Save_A3R(const char* file_name, std::vector<Atom2D>* as, LONG n, FLT r)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::out)) return FALSE;
#else
	if(!file.open(file_name, ios::out | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	header.count = n;
	header.r = r;
    std::vector<Atom2D>::iterator i;
	file.sputn((char*) &header, sizeof(header));
	for (i = as->begin(); i != as->end(); i++)
	{
		Vect3D r(i->r.x, i->r.y, 0);
	//	r.RotateX(DEG * .001);
	//	r.RotateY(DEG * .001);
		Vect3D::SFLOAT sr = r;
		file.sputn((char*) &sr, sizeof(Vect3D::SFLOAT));
	}

	file.close();

	return TRUE;
}

BOOL Save_A3R(const char* file_name, std::vector<Atom2D>* as, LONG n, FLT r, double w, double h)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::out)) return FALSE;
#else
	if(!file.open(file_name, ios::out | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	header.count = n;
	header.r = r;
    std::vector<Atom2D>::iterator i;
	file.sputn((char*) &header, sizeof(header));
	for (i = as->begin(); i != as->end(); i++)
	{
		Vect3D r(i->r.x-w*0.5, i->r.y-h*0.5, 0);
	//	r.RotateX(DEG * .001);
	//	r.RotateY(DEG * .001);
		Vect3D::SFLOAT sr = r;
		file.sputn((char*) &sr, sizeof(Vect3D::SFLOAT));
	}

	file.close();

	return TRUE;
}



//-------------------------------------------------------------------------------------------

BOOL SaveMarked_A3R(const char* file_name, Atom2D* as, LONG n, FLT r)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::out)) return FALSE;
#else
	if(!file.open(file_name, ios::out | ios::binary)) return FALSE;
#endif

	Atom2D* i;
	LONG n_marked = 0;
	for (i = as; i != as + n; i++) if(i->c) n_marked++;

	A3R_HEADER header;
	header.count = n_marked;
	header.r = r;

	file.sputn((char*) &header, sizeof(header));
	for (i = as; i != as + n; i++)
	{
		if(i->c)
		{
			Vect3D r(i->r.x, i->r.y, 0);
//			r.RotateX(DEG * .001);
//			r.RotateY(DEG * .001);
			Vect3D::SFLOAT sr = r;
			file.sputn((char*) &sr, sizeof(Vect3D::SFLOAT));
		}
	}

	file.close();

	return TRUE;
}

//-------------------------------------------------------------------------------------------

BOOL Get_count_A3R(const char* file_name, LONG& n)
{
	filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::in | ios::nocreate)) return FALSE;
#else
//	if(!file.open(file_name, ios::in | ios::nocreate | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	file.sgetn((char*)&header, sizeof(A3R_HEADER));
	if(strcmp(header.file_type, "a3r")) return FALSE;
	n = header.count;

	return TRUE;
}

//-------------------------------------------------------------------------------------------

BOOL Load_A3R(const char* file_name, Atom2D* as, LONG& n, FLT& r)
{
	/*filebuf file;

#ifdef UNIX
	if(!file.open(file_name, ios::in | ios::nocreate)) return FALSE;
#else
//	if(!file.open(file_name, ios::in | ios::nocreate | ios::binary)) return FALSE;
#endif

	A3R_HEADER header;
	file.sgetn((char*)&header, sizeof(A3R_HEADER));
	if(strcmp(header.file_type, "a3r")) return FALSE;
	n = min(n, header.count);
	r = header.r;

	for (Atom2D* i = as; i != as + n; i++)
	{
		Vect3D::SFLOAT sr;
		file.sgetn((char*) &sr, sizeof(Vect3D::SFLOAT));
		i->r.Set(sr.x, sr.y);
		i->v.Set(0, 0);
	}

	file.close();*/

	return TRUE;
}



//-------------------------------------------------------------------------------------------

