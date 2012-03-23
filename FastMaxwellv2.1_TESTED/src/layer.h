/***************************************************************************
 *   Copyright (C) 2005 by Xin Hu   *
 *   xinhu@mit.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef LAYER_H
#define LAYER_H

#include <complex>
#include <iostream>


class Layer
{
public:
	Layer(char *_name, double _epir, double _sigma, double _startH);
	
	char *name(void)const { return name_;};
	double startHeight(void)const {return startHeight_;};
	double epir(void) const {return epir_;};
	double sigma(void) const {return sigma_;};
	std::complex<double> k(void) const {return k_;};
	
	void set_epi(std::complex<double> epi) {epi_=epi;};
	void set_k(std::complex<double> k) {k_=k;};
private:
	char *name_;
	double epir_, sigma_, startHeight_;
	std::complex<double> epi_, k_;	
};

#endif
