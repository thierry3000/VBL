// compile with 
//
// g++-mp-6 -fopenmp -frounding-math -O2 -DNDEBUG -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lCGAL -lgmpxx -lmpfr -lgmp -lmpfi -lboost_thread-mt -lm test.cpp -o test
// g++-mp-7 -fopenmp -frounding-math -O2 -DNDEBUG -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lCGAL -lgmpxx -lmpfr -lgmp -lmpfi -lboost_thread-mt -lm test.cpp -o test


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>
#include <random>
#include <iostream>

	typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
	typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
	typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
	
	//Use the Fast_location tag. Default or Compact_location works too.
	
	typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
	typedef Delaunay::Point                                             Point;
	
	std::default_random_engine 											generator;
	std::uniform_real_distribution<double>	distribution(0.0,1.0);
	
int main()
{
	  int np=0;
	  
	  std::cout << "Enter number of points: ";
	  std::cin >> np;
	  
	  int nr=0;
	  
	  std::cout << "Enter number of repeats: ";
	  std::cin >> nr;
	  
	  for(int l=0; l < nr; l++)
	  	{
	  
		  std::vector< std::pair<Point,unsigned> > points;
	  
		  for(int k=0; k < np; k++)
			{
				double xc = distribution(generator);
				double yc = distribution(generator);
				double zc = distribution(generator);
				points.push_back( std::make_pair(Point(xc,yc,zc),k) );
				if( np < 10 ) std::cout << "(x,y,z) = (" << xc << "," << yc << "," << zc << ")" << std::endl;
				}
		
  
			  Delaunay T( points.begin(),points.end() );
	  
			  // CGAL_assertion( T.number_of_vertices() == 6 );
	  
			  // check that the info was correctly set.
			  Delaunay::Finite_vertices_iterator vit;
			  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
				if( points[ vit->info() ].first != vit->point() ){
				  std::cerr << "Error different info" << std::endl;
				  exit(EXIT_FAILURE);
				}
			  std::cout << l << "\t OK \r";
			  std::cout.flush();
	  
	  	}
	  
	  return 0;
}
