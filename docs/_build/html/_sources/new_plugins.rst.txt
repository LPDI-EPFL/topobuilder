.. _new_plugins:

Developing Your Own New Plugins
===============================

Although the TopoBuilder already contains a lot of useful Plugins, you may wish to create your own
module with wrapping around other external tools.

Implementing a plugin by yourself is straightforward and you may want to take a look at the code of
one of the base plugins to find out how the various plugins are implemented. However, we give a brief
introduction on what is needed to setup your own plugin.

Each plugin must be wrapped into a :class:`.Node` and two fundamental methods must be implemented:

* ``single_check()``:
  This function will do two things. A sanity check on a dummy example input that checks if all
  required information needed was passed, and instantiates the objects to save the calculated data to.

* ``single_execute()``:
  This function will actually do the calculation for a single :class:`.Case`. Note that multi_execute()
  has not been implemented yet.

Let's briefly take a look a an example from the make_topology plugin:

For the ``single_check()``, we first create a :class:`.Case` from saved dummy data. Then we check if
the required fields are present within the :class:`.Case` data. At the end we create an empty dictionary
where we will save our calculated data to to pass it on to the next plugin.

.. code-block:: python
  :linenos:

  def single_check( self, dummy: Dict ) -> Dict:
    kase = Case(dummy)

    # Check what it needs
    for itag in self.REQUIRED_FIELDS:
        if kase[itag] is None:
            raise NodeDataError(f'Field "{itag}" is required')

    # Include keywords
    kase.data.setdefault('metadata', {}).setdefault('equivalent_connectivities', [])
    return kase.data


For the ``single_execute()``, we load the actually current data into a :class:`.Case` and do
our computation, save the data and return it to pass on to the next plugin.

.. code-block:: python
  :linenos:
  
  def single_execute( self, data: Dict ) -> Dict:
      kase = Case(data)

      new_cases = []
      # If connectivities are pre-specified, only make those.
      if kase.connectivity_count > 0:
          new_cases.extend(self.eval_representatives(kase, self.representatives, self.sampling))
      else:
          new_cases.extend(self.eval_representatives(
                           self.explore_connectivities(kase), self.representatives, self.sampling))
      self.log.notice(f'case count: {len(new_cases)}')
      return new_cases[0]
