#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  json.py
"""
Custom json encoder to support PyMassSpec NIST Search classes
"""
#
#  This file is part of PyMassSpec NIST Search
#  Python interface to the NIST MS Search DLL
#
#  Copyright (c) 2020 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  PyMassSpec NIST Search is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#
#  PyMassSpec NIST Search is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  PyMassSpec NIST Search includes the redistributable binaries for NIST MS Search in
#  the x86 and x64 directories. Available from
#  ftp://chemdata.nist.gov/mass-spc/v1_7/NISTDLL3.zip .
#  ctnt66.dll and ctnt66_64.dll copyright 1984-1996 FairCom Corporation.
#  "FairCom" and "c-tree Plus" are trademarks of FairCom Corporation
#  and are registered in the United States and other countries.
#  All Rights Reserved.


# 3rd party
from pyms.json import PyMassSpecEncoder
from sdjson import register_encoder

# this package
from .base import NISTBase
from .reference_data import ReferenceData
from .search_result import SearchResult


class PyNISTEncoder(PyMassSpecEncoder):
    """
    Custom json encoder to support PyMassSpec NIST Search classes
    """

    def default(self, o):
        if isinstance(o, ReferenceData):
            return o.__dict__(recursive=True)
        elif isinstance(o, NISTBase):
            return dict(o)
        else:
            return PyMassSpecEncoder.default(self, o)


@register_encoder(ReferenceData)
def encode_reference_data(obj):
    return dict(obj)


@register_encoder(SearchResult)
def encode_search_result(obj):
    return dict(obj)
