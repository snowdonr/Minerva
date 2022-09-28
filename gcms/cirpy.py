# -*- coding: utf-8 -*-
"""
CIRpy
Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.
https://github.com/mcs07/CIRpy
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import functools
import inspect
import logging
import os

try:
    from urllib.error import HTTPError
    from urllib.parse import quote, urlencode
    from urllib.request import urlopen
except ImportError:
    from urllib import urlencode
    from urllib2 import quote, urlopen, HTTPError

try:
    from lxml import etree
except ImportError:
    try:
        import xml.etree.cElementTree as etree
    except ImportError:
        import xml.etree.ElementTree as etree


__author__ = 'Matt Swain'
__email__ = 'm.swain@me.com'
__version__ = '1.0.2'
__license__ = 'MIT'

log = logging.getLogger('cirpy')
log.addHandler(logging.NullHandler())

API_BASE = 'https://cactus.nci.nih.gov/chemical/structure'
FILE_FORMATS = {
    'alc', 'cdxml', 'cerius', 'charmm', 'cif', 'cml', 'ctx', 'gjf', 'gromacs', 'hyperchem', 'jme', 'maestro', 'mol',
    'mol2', 'mrv', 'pdb', 'sdf3000', 'sln', 'xyz'
}


def construct_api_url(compound_input, representation, resolvers=None, get3d=False, tautomers=False, xml=True, **kwargs):
    """Return the URL for the desired API endpoint.
    :param string compound_input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(str) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :param bool xml: (Optional) Whether to return full XML response
    :returns: CIR API URL
    :rtype: str
    """
    # File formats require representation=file and the format in the querystring
    if representation in FILE_FORMATS:
        kwargs['format'] = representation
        representation = 'file'
    # Prepend compound_input with 'tautomers:' to return all tautomers
    if tautomers:
        compound_input = 'tautomers:%s' % compound_input
    url = '%s/%s/%s' % (API_BASE, quote(compound_input), representation)
    if xml:
        url += '/xml'
    if resolvers:
        kwargs['resolver'] = ','.join(resolvers)
    if get3d:
        kwargs['get3d'] = True
    if kwargs:
        url += '?%s' % urlencode(kwargs)
    return url


def request(compound_input, representation, resolvers=None, get3d=False, tautomers=False, **kwargs):
    """Make a request to CIR and return the XML response.
    :param string compound_input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :returns: XML response from CIR
    :rtype: Element
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    url = construct_api_url(compound_input, representation, resolvers, get3d, tautomers, **kwargs)
    log.debug('Making request: %s', url)
    response = urlopen(url)
    return etree.parse(response).getroot()


class Result(object):
    """A single result returned by CIR."""

    def __init__(self, compound_input, notation, input_format, resolver, representation, value):
        """
        :param string compound_input: Originally supplied input identifier that produced this result
        :param string notation: Identifier matched by the resolver or tautomer ID
        :param string input_format: Format of the input as interpreted by the resolver
        :param string resolver: Resolver used to produce this result
        :param string representation: Requested output representation
        :param value: Actual result value
        :type value: string or list(string)
        """
        self.compound_input = compound_input
        self.representation = representation
        self.resolver = resolver
        self.input_format = input_format
        self.notation = notation
        self.value = value

    def __repr__(self):
        return 'Result(input=%r, representation=%r, resolver=%r, input_format=%r, notation=%r, value=%r)' \
               % (self.compound_input, self.representation, self.resolver, self.input_format, self.notation, self.value)

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.__dict__ == other.__dict__

    def __getitem__(self, prop):
        """Allow dict-style access to attributes to ease transition from when results were dicts."""
        if prop in self.__dict__:
            return getattr(self, prop)
        raise KeyError(prop)

    def __setitem__(self, prop, val):
        """Allow dict-style setting of attributes to ease transition from when results were dicts."""
        setattr(self, prop, val)

    def __contains__(self, prop):
        """Allow dict-style checking of attributes to ease transition from when results were dicts."""
        return prop in self.__dict__

    def to_dict(self):
        """Return a dictionary containing Result data."""
        return self.__dict__


def query(compound_input, representation, resolvers=None, get3d=False, tautomers=False, **kwargs):
    """Get all results for resolving input to the specified output representation.
    :param string compound_input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :returns: List of resolved results
    :rtype: list(Result)
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    tree = request(compound_input, representation, resolvers, get3d, tautomers, **kwargs)
    results = []
    for data in tree.findall('.//data'):
        value = [item.text for item in data.findall('item')]
        result = Result(
            compound_input=tree.attrib['string'],
            representation=tree.attrib['representation'],
            resolver=data.attrib['resolver'],
            input_format=data.attrib['string_class'],
            notation=data.attrib['notation'],
            value=value[0] if len(value) == 1 else value
        )
        results.append(result)
    log.debug('Received %s query results', len(results))
    return results


def resolve(compound_input, representation, resolvers=None, get3d=False, **kwargs):
    """Resolve input to the specified output representation.
    :param string compound_input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :returns: Output representation or None
    :rtype: string or None
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    # Take first result from XML query
    results = query(compound_input, representation, resolvers, False, get3d, **kwargs)
    result = results[0].value if results else None
    return result


def resolve_image(compound_input, resolvers=None, fmt='png', width=300, height=300, frame=False, crop=None, bgcolor=None,
                  atomcolor=None, hcolor=None, bondcolor=None, framecolor=None, symbolfontsize=11, linewidth=2,
                  hsymbol='special', csymbol='special', stereolabels=False, stereowedges=True, header=None, footer=None,
                  **kwargs):
    """Resolve input to a 2D image depiction.
    :param string compound_input: Chemical identifier to resolve
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param string fmt: (Optional) gif or png image format (default png)
    :param int width: (Optional) Image width in pixels (default 300)
    :param int height: (Optional) Image height in pixels (default 300)
    :param bool frame: (Optional) Whether to show border frame (default False)
    :param int crop: (Optional) Crop image with specified padding
    :param int symbolfontsize: (Optional) Atom label font size (default 11)
    :param int linewidth: (Optional) Bond line width (default 2)
    :param string bgcolor: (Optional) Background color
    :param string atomcolor: (Optional) Atom label color
    :param string hcolor: (Optional) Hydrogen atom label color
    :param string bondcolor: (Optional) Bond color
    :param string framecolor: (Optional) Border frame color
    :param bool hsymbol: (Optional) Hydrogens: all, special or none (default special)
    :param bool csymbol: (Optional) Carbons: all, special or none (default special)
    :param bool stereolabels: (Optional) Whether to show stereochemistry labels (default False)
    :param bool stereowedges: (Optional) Whether to show wedge/dash bonds (default True)
    :param string header: (Optional) Header text above structure
    :param string footer: (Optional) Footer text below structure
    """
    # Aggregate all arguments into kwargs
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    for arg in args:
        if values[arg] is not None:
            kwargs[arg] = values[arg]
    # Turn off anti-aliasing for transparent background
    if kwargs.get('bgcolor') == 'transparent':
        kwargs['antialiasing'] = False
    # Renamed parameters
    if 'stereolabels' in kwargs:
        kwargs['showstereo'] = kwargs.pop('stereolabels')
    if 'fmt' in kwargs:
        kwargs['format'] = kwargs.pop('fmt')
    # Toggle stereo wedges
    if 'stereowedges' in kwargs:
        status = kwargs.pop('stereowedges')
        kwargs.update({'wedges': status, 'dashes': status})
    # Constant values
    kwargs.update({'representation': 'image', 'xml': False})
    url = construct_api_url(**kwargs)
    log.debug('Making image request: %s', url)
    response = urlopen(url)
    return response.read()


# TODO: Support twirl as fmt paramter?
# TODO: ipython html repr twirl, ipython png repr image


def download(compound_input, filename, representation, overwrite=False, resolvers=None, get3d=False, **kwargs):
    """Convenience function to save a CIR response as a file.
    This is just a simple wrapper around the resolve function.
    :param string compound_input: Chemical identifier to resolve
    :param string filename: File path to save to
    :param string representation: Desired output representation
    :param bool overwrite: (Optional) Whether to allow overwriting of an existing file
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    :raises IOError: if overwrite is False and file already exists
    """
    result = resolve(compound_input, representation, resolvers, get3d, **kwargs)
    # Just log and return if nothing resolved
    if not result:
        log.debug('No file to download.')
        return
    # Only overwrite an existing file if explicitly instructed to.
    if not overwrite and os.path.isfile(filename):
        raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % filename)
    # Ensure file ends with a newline
    if not result.endswith('\n'):
        result += '\n'
    with open(filename, 'w') as f:
        f.write(result)


def memoized_property(fget):
    """Decorator to create memoized properties."""
    attr_name = '_{0}'.format(fget.__name__)

    @functools.wraps(fget)
    def fget_memoized(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fget(self))
        return getattr(self, attr_name)
    return property(fget_memoized)


class Molecule(object):
    """Class to hold and cache the structure information for a given CIR input."""

    def __init__(self, compound_input, resolvers=None, get3d=False, **kwargs):
        """Initialize with a resolver input."""
        self.compound_input = compound_input
        self.resolvers = resolvers
        self.get3d = get3d
        self.kwargs = kwargs
        log.debug('Instantiated Molecule: %s' % self)

    def __repr__(self):
        return 'Molecule(input=%r, resolvers=%r, get3d=%r, kwargs=%r)' \
               % (self.compound_input, self.resolvers, self.get3d, self.kwargs)

    @memoized_property
    def stdinchi(self):
        """Standard InChI."""
        return resolve(self.compound_input, 'stdinchi', self.resolvers, **self.kwargs)

    @memoized_property
    def stdinchikey(self):
        """Standard InChIKey."""
        return resolve(self.compound_input, 'stdinchikey', self.resolvers, **self.kwargs)

    @memoized_property
    def inchi(self):
        """Non-standard InChI. (Uses options DONOTADDH W0 FIXEDH RECMET NEWPS SPXYZ SAsXYZ Fb Fnud)."""
        return resolve(self.compound_input, 'inchi', self.resolvers, **self.kwargs)

    @memoized_property
    def smiles(self):
        """SMILES string."""
        return resolve(self.compound_input, 'smiles', self.resolvers, **self.kwargs)

    @memoized_property
    def ficts(self):
        """FICTS NCI/CADD hashed structure identifier."""
        return resolve(self.compound_input, 'ficts', self.resolvers, **self.kwargs)

    @memoized_property
    def ficus(self):
        """FICuS NCI/CADD hashed structure identifier."""
        return resolve(self.compound_input, 'ficus', self.resolvers, **self.kwargs)

    @memoized_property
    def uuuuu(self):
        """uuuuu NCI/CADD hashed structure identifier."""
        return resolve(self.compound_input, 'uuuuu', self.resolvers, **self.kwargs)

    @memoized_property
    def hashisy(self):
        """CACTVS HASHISY identifier."""
        return resolve(self.compound_input, 'hashisy', self.resolvers, **self.kwargs)

    @memoized_property
    def sdf(self):
        """SDF file."""
        return resolve(self.compound_input, 'sdf', self.resolvers, **self.kwargs)

    @memoized_property
    def names(self):
        """List of chemical names."""
        return resolve(self.compound_input, 'names', self.resolvers, **self.kwargs)

    @memoized_property
    def iupac_name(self):
        """IUPAC approved name."""
        return resolve(self.compound_input, 'iupac_name', self.resolvers, **self.kwargs)

    @memoized_property
    def cas(self):
        """CAS registry numbers."""
        return resolve(self.compound_input, 'cas', self.resolvers, **self.kwargs)

    @memoized_property
    def mw(self):
        """Molecular weight."""
        return resolve(self.compound_input, 'mw', self.resolvers, **self.kwargs)

    @memoized_property
    def formula(self):
        """Molecular formula"""
        return resolve(self.compound_input, 'formula', self.resolvers, **self.kwargs)

    @memoized_property
    def h_bond_donor_count(self):
        """Hydrogen bond donor count."""
        return resolve(self.compound_input, 'h_bond_donor_count', self.resolvers, **self.kwargs)

    @memoized_property
    def h_bond_acceptor_count(self):
        """Hydrogen bond acceptor count."""
        return resolve(self.compound_input, 'h_bond_acceptor_count', self.resolvers, **self.kwargs)

    @memoized_property
    def h_bond_center_count(self):
        """Hydrogen bond center count."""
        return resolve(self.compound_input, 'h_bond_center_count', self.resolvers, **self.kwargs)

    @memoized_property
    def rule_of_5_violation_count(self):
        """Rule of 5 violation count."""
        return resolve(self.compound_input, 'rule_of_5_violation_count', self.resolvers, **self.kwargs)

    @memoized_property
    def rotor_count(self):
        """Rotor count."""
        return resolve(self.compound_input, 'rotor_count', self.resolvers, **self.kwargs)

    @memoized_property
    def effective_rotor_count(self):
        """Effective rotor count."""
        return resolve(self.compound_input, 'effective_rotor_count', self.resolvers, **self.kwargs)

    @memoized_property
    def ring_count(self):
        """Ring count."""
        return resolve(self.compound_input, 'ring_count', self.resolvers, **self.kwargs)

    @memoized_property
    def ringsys_count(self):
        """Ring system count."""
        return resolve(self.compound_input, 'ringsys_count', self.resolvers, **self.kwargs)

    @memoized_property
    def image(self):
        """2D image depiction."""
        return resolve_image(self.compound_input, self.resolvers, **self.kwargs)

    @property
    def image_url(self):
        """URL of a GIF image."""
        return construct_api_url(self.compound_input, 'image', self.resolvers, False, self.get3d, False, **self.kwargs)

    @property
    def twirl_url(self):
        """Url of a TwirlyMol 3D viewer."""
        return construct_api_url(self.compound_input, 'twirl', self.resolvers, False, self.get3d, False, **self.kwargs)

    def download(self, filename, representation, overwrite=False):
        """Download the resolved structure as a file.
        :param string filename: File path to save to
        :param string representation: Desired output representation
        :param bool overwrite: (Optional) Whether to allow overwriting of an existing file
        """
        download(self.compound_input, filename, representation, overwrite, self.resolvers, self.get3d, **self.kwargs)
