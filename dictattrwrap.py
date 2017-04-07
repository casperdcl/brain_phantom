class DictAttrWrap(object):
  def __init__(self, *a, **k):
    self.d = dict(*a, **k)

  def __getattr__(self, k):
    return self.d.get('--' + k, self.d.get('<' + k + '>'))
