import numpy as N
import new, functools
from functools import partial

"""
Stuff not directly related to the FEM code. Stolen from all over, probably, and
not well tested by me, directly at least.
"""
class memoized(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.

    From http://wiki.python.org/moin/PythonDecoratorLibrary#head-11870a08b0fa59a8622201abfac735ea47ffade5
    """
    def __init__(self, func):
       self.func = func
       self.cache = {}
    def __call__(self, *args):
       try:
          return self.cache[args]
       except KeyError:
          self.cache[args] = value = self.func(*args)
          return value
       except TypeError:
          # uncachable -- for instance, passing a list as an argument.
          # Better to not cache than to blow up entirely.
          return self.func(*args)

def memoized2(func):
    cache = {}
    def cached_fun(*args):
       try:
          return cache[args]
       except KeyError:
          cache[args] = value = func(*args)
          return value
       except TypeError:
          # uncachable -- for instance, passing a list as an argument.
          # Better to not cache than to blow up entirely.
          return func(*args)

    return cached_fun

class memoize_last(object):
    """
    Decorator that caches the last function call.
    """

    def __init__(self,func):
        self.func = func
        self.last_args = None

    def __call__(self, *args):
        if self.last_args == args:
            return self.cached
        else:
            self.cached = value = self.func(*args)
            self.last_args = args
            return value

    def clearCache(self):
        self.last_args = None
        try: del(self.cached)
        except AttributeError: pass

    def __getattr__(self, attr):
        return getattr(self.func, attr)

class CacheLast(type):
    """
    Alternative to method_memoize_last that applies the standard memoize_last
    decorator at object instantiation time. Set __metaclass__ = CacheLast and
    decorate functions that are to be cached with @CacheLast.CachedMethod
    """
    def __call__(cls, *names, **kwargs):
        obj = super(CacheLast, cls).__call__(*names, **kwargs)
        obj.cachedAttrs = dict()
        # Loop over all the class attributes
        for attrname, attrobj in filter_attrs(
            cls, lambda attrname, attrobj: isinstance(attrobj, cls.CachedMethod)):
            # Get the unbound method out of the CachedMethod object and
            # bind it to the object (i.e. instance)
            cachedmethod = memoize_last(
                new.instancemethod(attrobj.meth, obj, cls))
            setattr(obj, attrname, cachedmethod)
            obj.cachedAttrs[attrname] = cachedmethod
        return obj

    class CachedMethod(object):
        def __init__(self, meth):
            self.meth = meth


class method_memoize_last(object):
    """
    Decorator that caches the last method call.
    
    Based in part on the http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/325205
    comment:
    
    A working example, Chris Spencer, 2006/04/27
    This will corretly decorate a method:
`    """
    def __init__(self, func):
        self.func = func
        self.instances={}
    
    def __get__(self, instance, cls=None):
        self.instance = instance
        return self

    def __call__(self,*args):
        instance = self.instance
        instances = self.instances
        try:
            (last_args, cached) = instances[instance]
        # This instance has not been cached yet
        except KeyError:
            value = self.func(instance, *args)
            instances[instance]=(args, value)
            return value
        # This instance has been cached, compare method arguments
        if last_args == args:
            return cached
        else:
            value = self.func(instance, *args)
            instances[instance] = (args, value)
            return value

def old_partial(*args, **kwargs):
    """
    From comments of
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52549 I renamed it
    to partial from curry, since its use seems identical to the partial
    function proposed for Python 2.5

    use:

    fp = partial(func, parm1, kwarg=val)

    fp(x) is now equivalent to f(parm1, x, kwarg=val)
    
    """
    function, args = args[0], args[1:]
    if args and kwargs:
        def result(*rest, **kwrest):
            combined = kwargs.copy()
            combined.update(kwrest)
            return function(*args + rest, **combined)
    elif args:
        if len(args) > 1 or args[0] is None:
            def result(*rest, **kwrest):
                return function(*args + rest, **kwrest)
        else:
            # Special magic: make a bound object method on the arg
            return new.instancemethod(function, args[0], object)
    elif kwargs:
        def result(*rest, **kwrest):
            if kwrest:
                combined = kwargs.copy()
                combined.update(kwrest)
            else:
                combined = kwargs
            return function(*rest, **combined)
    else:
        return function
    return result

class Struct(dict):
    forbidden = set(dir(dict) + ['forbidden'])
    def __setitem__(self, key, item):
        try:
            self.__setattr__(key, item)
        except AttributeError:
            raise KeyError,  "Illegal key"

    def __delattr__(self, name):
        del self[name]
        dict.__delattr__(self, name)

    def __init__(self, *names, **kwargs):
        dict.__init__(self, *names, **kwargs)
        for k,v in self.iteritems():
            assert(k not in self.forbidden)
            dict.__setattr__(self, k, v)

    def update(self, in_dict):
        for k,v in in_dict.iteritems():
            setattr(self, k, v)

    def __setattr__(self, name, value):
        if name in self.forbidden: raise AttributeError, "Illegal atribute name"
        dict.__setattr__(self, name, value)
        dict.__setitem__(self, name, value)
    

isBound = lambda meth: meth.im_self is not None

def add_attr(attr, val):
    def add_attr_(fun):
        setattr(fun, attr, val)
        return fun
    return add_attr_

def filter_attrs(obj, cond):
    return ((attrname, attrobj) for attrname, attrobj in
            map(lambda key: (key, getattr(obj, key)), dir(obj))
            if cond(attrname, attrobj))

def max_or_scalar(x):
    try: return N.max(x)
    except TypeError: return x

from itertools import chain
class rechain(object):
    def __init__(self, *names):
        self.names = names
    def __iter__(self):
        return chain(*self.names)

def close_to_point(point, eps):
    def close_to_pointp(x):
        return N.abs(point - x) < eps
    return close_to_pointp

def almost_leq(eps):
    def _leq(lesser, greater):
        return lesser - greater < eps
    return _leq

def almost_geq(eps):
    def _geq(greater, lesser):
        return greater - lesser > -eps
    return _geq

def in_box(bound_n, bound_p):
    def in_boxp(x):
        return (N.all(x <= bound_p) and N.all(x >= bound_n))
    return in_boxp

def in_box_vec(bound_n, bound_p):
    x_n, y_n, z_n = bound_n
    x_p, y_p, z_p = bound_p
    def in_box_vecp(r):
        x,y,z = r.T
        return (((x <= x_p) & (x >= x_n)) &
                ((y <= y_p) & (y >= y_n)) &
                ((z <= z_p) & (z >= z_n)))
    return in_box_vecp

def on_box_surf(bound_n, bound_p, eps):
    xp_p, yp_p, zp_p = [close_to_point(u, eps) for u in bound_p]
    xn_p, yn_p, zn_p = [close_to_point(u, eps) for u in bound_n]
    xp, yp, zp = bound_p
    xn, yn, zn = bound_n
    def on_box_surfp(r):
        x,y,z = r.T
        return ( (N.all(xn_p(x)) or N.all(xp_p(x))) and N.all(
            (yn <= y) & (y <= yp) & (zn <= z) & (z <= zp)) or
                 (N.all(yn_p(y)) or N.all(yp_p(y))) and N.all(
            (zn <= z) & (z <= zp) & (xn <= x) & (x <= xp)) or
                 (N.all(zn_p(z)) or N.all(zp_p(z))) and N.all(
            (xn <= x) & (x <= xp) & (yn <= y) & (y <= yp)) )
                 
    return on_box_surfp

class ArrayMemory(object):
    def __init__(self, model_array, mem_len=3):
        self.mem = [N.zeros_like(model_array) for i in range(mem_len)]

    def push(self, val):
        mem = self.mem
        assert(val.shape == mem[0].shape)
        assert(val.dtype is mem[0].dtype)
        mem.insert(0, val.copy())
        mem.pop()

    def __getitem__(self, index):
        return self.mem[index]
    
    def __iter__(self):
        return self.mem.__iter__()

def NoArgsCachedMethod(meth):
    cacheAttr = '_'+meth.func_name
    def c_meth(self):
        try: return getattr(self, cacheAttr)
        except AttributeError:
            val = f(self)
            setattr(self, cacheAttr, val)
        return val
    return functools.update_wrapper(c_meth, meth)

class attr_poly1d(N.poly1d):
    def __setattr__(self, key, val):
        self.__dict__[key] = val

def gen_lagrange_polys(points):
    def make_poly(int_pt, zero_pts):
        poly = attr_poly1d(N.poly(zero_pts)/N.multiply.reduce(
            [int_pt - p for p in zero_pts]))
        poly.pt = int_pt
        return poly
    return [make_poly(pi, [pz for pz in points if pz != pi]) for pi in points]

RMS = lambda x: N.sqrt(N.sum(x**2)/len(x))
RMS_weighted = lambda x, w: N.sqrt(N.sum((x**2)*w)/N.sum(w))

def find_peaks(vals, cutoff=0):
    c_vals = vals[1:-1]
    l_vals = vals[0:-2]
    r_vals = vals[2:]
    peaks = (c_vals > l_vals) & (c_vals > r_vals) & ( c_vals >= cutoff)
    peak_ind = N.arange(len(c_vals))[peaks] +1
    return peak_ind

def remove_neg(arr, copy=True):
    return N.array(arr[arr >= 0], copy=copy)

def remove_constr(entnos, entlist, copy=True):
    return N.array(entnos[entlist[entnos].isFree], copy=copy)

class ProxyArray(object):
    __getitem__ = lambda self, i: self.subarr[i]
    __iter__ = lambda self: self.subarr.__iter__()
    __array__ = lambda self: self.subarr.__array__()
    __len__ = lambda self: self.subarr.__len__()

