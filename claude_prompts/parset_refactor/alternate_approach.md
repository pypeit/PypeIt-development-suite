# Assessment: `ParSet`/`pypeitpar.py` design vs. prior art, and proposals for two corner cases

## Context

The `parset_refactor` branch replaced the old imperative `ParSet` pattern
(parallel `defaults`/`options`/`dtypes`/`descr` dicts matched by list
position) with a declarative one: each subclass defines a class-level
`parameters = {'key': parset.set_parameter_definition(dtype=..., default=...,
options=..., descr=...), ...}` dict, and the generic `ParSet` base class
handles construction, dict-style access, dtype checking, and all I/O
(`to_dict`/`from_dict`/`to_config`/`to_header`/`from_header`). This is a
substantial, hand-rolled piece of infrastructure, and the user does not want
to redo it again soon — so before building further on top of it (the
`funcpar_update`/`sensfunc_merge` branches), it's worth checking the design
against established patterns elsewhere in the Python ecosystem and hardening
two corner cases that came up organically during the refactor:

1. **Callable defaults** — `ReduxPar.redux_path` needs a default that can
   only be computed at instantiation time (`Path.cwd()`), not at class
   definition/import time. The current fix: include `Callable` in `dtype`,
   store the unevaluated callable as the value, and rely on
   `ParSet.fill_callable()` (invoked from `PypeItPar.validate()`) to replace
   it with the call result later.
2. **Signature-derived parameter sets** (`pypeit/par/funcpar.py::FuncPar`,
   prototyped on `funcpar_update` and extended on `sensfunc_merge` with
   `DifferentialEvolutionPar`/`DJSRejectPar`) — a `ParSet` subclass whose
   `parameters` dict is built automatically from a wrapped function's
   keyword-argument signature (via a custom metaclass,
   `FuncParMetaClass`), so the exposed parameters can't drift from the
   real function API. The user's own comment in that file already flags
   discomfort with needing a metaclass for this.

## Current mechanics (as implemented on this branch)

- `set_parameter_definition(dtype=None, default=None, options=None,
  descr=None)` (`pypeit/par/parset.py:56-111`) just builds a plain dict;
  there's no separate "how to compute the default" concept — `default` is
  always a value, and "is this a callable-producing-a-default" is inferred
  purely by including `collections.abc.Callable` in `dtype`.
- `ParSet.__init__` (`parset.py:167-214`) sets every key to
  `self.parameters[key]['default']` through the normal `__setitem__` dtype
  check, then applies user kwargs, then calls `self.validate()`.
  `__setitem__` (`parset.py:227-311`) does `isinstance(value, d)` against
  `dtype`; since `Callable` is an ABC matching anything with `__call__`,
  an unevaluated `Path.cwd` (a bound method) passes and is stored *as the
  callable itself*.
- Resolution only happens via `ParSet.fill_callable()`
  (`parset.py:550-574`), which recurses into nested `ParSet`s and replaces
  any `callable(self[key])` value with `self[key]()`. This is invoked from
  `PypeItPar.validate()` (`pypeitpar.py:5373-5375`) — i.e. **only** when a
  full `PypeItPar` tree is constructed. A standalone `ReduxPar()` (or any
  nested `ParSet` built and used on its own, e.g. in a unit test) leaves
  `redux_path` as a raw bound method until/unless something explicitly
  calls `fill_callable()`. Since `to_dict`/`to_config`/`to_header`
  (`parset.py:691-1003`) just serialize whatever is in `_data`, this is a
  latent bug: any I/O performed on such an instance before `fill_callable()`
  runs would emit a method's `repr()`/`str()` instead of a real path.
- `FuncPar` (`funcpar_update`/`sensfunc_merge` branches,
  `pypeit/par/funcpar.py`) uses `FuncParMetaClass.__new__` to introspect
  `cls.func`'s keyword arguments (via `utils.get_func_kwargs`, presumably
  wrapping `inspect.signature`) at class-definition time, filters them
  through `kw_subset`/`omitted_keys`, and assigns the resulting
  `parameters` dict as a class attribute. `dtype` is left unconstrained
  except that `sensfunc_merge` forces `dtype=tuple` when a keyword's
  default is itself a tuple. Neither branch gives `FuncPar` a way to
  actually call the wrapped function with the collected values — that's
  left to the caller.
- Nested/dict-style access (`par['calibrations']['slitedges']['pad']`) and
  the `FrameGroupPar.parameters | {...}` merge idiom for subclassing both
  work well today and should be kept as-is.

## Comparison to prior art

| Pattern | Declarative fields | Lazy/computed defaults | Dict-style (`par['x']`) access | Signature-derived fields | Extra dependency |
|---|---|---|---|---|---|
| **This branch's `ParSet`** | Yes (`parameters` dict) | Ad hoc (`Callable` in `dtype` + manual `fill_callable()`) | Yes, native | Ad hoc (metaclass) | None |
| `dataclasses` (stdlib) | Yes (`field(...)`) | Yes — `field(default_factory=...)`, the direct precedent for corner case 1 | No (attribute access only) | No | None |
| `attrs`/`cattrs` | Yes (`attr.ib(...)`) | Yes — `attr.ib(factory=...)`, plus per-field validators/converters | No | No | Yes |
| `pydantic` v2 | Yes (`Field(...)`) | Yes — `Field(default_factory=...)` | No natively (attribute access; dict-style needs a mixin) | Yes — `create_model`/`validate_call` derive a model from a function signature, the closest precedent for corner case 2 | Yes (heavier) |
| `traitlets` (IPython/Jupyter config; already an indirect dependency via `IPython.embed`) | Yes (`Int()`, `Instance()`, `Callable()` trait types) | Yes — `@default('name')` dynamic-default decorator, resolved lazily on first *access*, not eagerly at construction | No | No | Already present transitively |
| `sklearn.BaseEstimator` | No fixed schema | N/A | No | Yes — `get_params()`/`set_params()` introspect `type(self).__init__`'s signature **at call time**, via a plain method, no metaclass at all | Not applicable (sklearn-specific idiom, not a library to depend on) |
| OmegaConf/Hydra structured configs | Yes (dataclass-backed) | Via custom resolvers (`${...}`) | Yes, native | No | Yes (heavy, YAML-centric) |

Takeaway: nothing here beats keeping the hand-rolled `ParSet` outright — it's
the only option with dict-style nested access *and* the FITS-header I/O
PypeIt needs, dependency-free. But two idioms are worth borrowing directly
into it rather than re-deriving ad hoc solutions:

- **`dataclasses.field(default_factory=...)`** is the standard, minimal
  answer to corner case 1.
- **`sklearn.BaseEstimator`'s metaclass-free signature introspection**
  and **Python's `__init_subclass__` hook (PEP 487, added in Python 3.6 —
  well below PypeIt's `requires-python = ">=3.11,<3.15"` floor, so no
  version concern)** together answer the user's own discomfort with
  `FuncParMetaClass` in corner case 2 — `__init_subclass__` is the
  language's sanctioned replacement for exactly this "customize subclass
  creation" use case, without the risks of a custom metaclass (harder to
  read, and it would conflict with any other metaclass `ParSet`/`FuncPar`
  might need later, e.g. `ABCMeta`).

## Proposal 1 — explicit `default_factory` for lazy/computed defaults

- Add a `default_factory` kwarg to `set_parameter_definition()`, mutually
  exclusive with `default` (raise if both are given). Drop `Callable` from
  `dtype` for this purpose entirely — `dtype` should describe only the
  *resolved* value's type.
- In `ParSet.__init__`, when setting a key's initial value, call
  `self.parameters[key]['default_factory']()` if present instead of using
  `default` directly. This resolves the value **eagerly, at construction
  time**, for every `ParSet`/subclass — not just ones reachable from a
  top-level `PypeItPar.validate()` cascade — closing the I/O-leak gap
  described above.
- Remove `fill_callable()` and its call from `PypeItPar.validate()`
  entirely; it becomes unnecessary since every nested `ParSet` resolves its
  own factory-based defaults in its own `__init__`.
- `ReduxPar.redux_path` becomes:
  `dtype=[str, Path], default_factory=Path.cwd, descr=...` — no more
  `Callable` in the type list, and the stored value is always a concrete
  `Path`/`str` from the moment the instance exists.
- **Advantage over current**: never possible to observe or serialize an
  unresolved callable; matches a well-known stdlib idiom, easy to document
  by reference to `dataclasses`. **Disadvantage**: purely eager — loses the
  theoretical ability to defer the computation until the value is actually
  read (not a real cost for `Path.cwd()`, but worth naming as the general
  trade-off versus e.g. `traitlets`' access-time-lazy `@default`).

## Proposal 2 — `FuncPar` without a metaclass, plus signature-derived `dtype`

- Replace `FuncParMetaClass.__new__` with a `FuncPar.__init_subclass__(cls,
  **kwargs)` classmethod: same effect (computing and attaching `parameters`
  when a subclass like `DifferentialEvolutionPar` is defined), no metaclass,
  no metaclass-conflict risk, and it directly answers "I'm not crazy about
  needing a metaclass" with the standard modern alternative.
- While touching `_define_parameters`, use `inspect.signature(func)` /
  `typing.get_type_hints(func)` to populate `dtype` from the wrapped
  function's own type annotations wherever one exists and is a plain
  concrete type (or a simple `Optional[T]`/`Union` collapsible to a
  `dtype` list) — generalizing the tuple-only heuristic already added on
  `sensfunc_merge`, so parameters like `strategy: str = 'best1bin'` get a
  real `dtype=[str]` instead of `None`. Fall back to `None` (unconstrained)
  when annotations are absent or too complex to translate.
- Out of scope: giving `FuncPar` a way to actually invoke the wrapped
  function (e.g. a `__call__` method) — the user wants `FuncPar` to remain
  limited to exposing/collecting the parameter values; calling `self.func`
  with them is left to the caller, as it already is on both
  `funcpar_update` and `sensfunc_merge`.
- **Alternative considered and not recommended**: mirror
  `sklearn.BaseEstimator` even more closely by making `parameters` a
  lazily-computed property (introspecting `cls.func` on first access)
  instead of a class attribute populated at subclass-creation time. This
  needs no `__init_subclass__` at all, but it breaks the convention, used
  everywhere else in `pypeitpar.py`, that `SomeClass.parameters` is a plain
  dict literal inspectable without instantiation (relied on by
  `doc/scripts/build_par_rst.py` and the `FrameGroupPar.parameters | {...}`
  merge idiom) — so `__init_subclass__` is the better fit here, keeping
  `parameters` a real, eagerly-populated dict.

## Advantages / disadvantages summary (relative to the current branch)

| | Current (this branch) | Proposed |
|---|---|---|
| Callable defaults | Ambiguous `dtype` (type-or-factory), deferred resolution only reachable via `PypeItPar.validate()`, silent I/O-leak risk for standalone nested `ParSet`s | Explicit `default_factory`, resolved eagerly and uniformly at every `ParSet.__init__`, no leak risk, matches `dataclasses` idiom |
| `FuncPar` construction | Custom metaclass (author-flagged discomfort), no signature-derived `dtype` beyond a tuple special case | `__init_subclass__` (standard, conflict-free, no minimum-version concern), annotation-derived `dtype`; still no invocation mechanism, unchanged from current scope |
| Risk of the change | — | Low: both changes are localized to `parset.py`/`funcpar.py` machinery; no change to `parameters = {...}` declarations elsewhere except `ReduxPar.redux_path`'s one entry, and no change to the dict-style access pattern the user wants preserved |

## Pydantic feasibility exploration (investigated, not adopted this pass)

`pydantic>=2.0` turns out to already be a project dependency, actively used
in `pypeit/state/*.py` and `pypeit/dashboard/model.py` (`BaseModel`,
`Field(default_factory=...)`, `Literal[...]` for restricted-value fields,
`model_dump_json`/`model_validate`) — so it is not a new dependency cost,
and it natively solves corner case 1 (`Field(default_factory=Path.cwd)`
needs no custom machinery at all) and fits corner case 2 well (a custom
`FuncParMetaClass` would additionally need to subclass pydantic's own
`ModelMetaclass` to avoid a metaclass conflict, whereas `__init_subclass__`
has no such conflict — another concrete reason to prefer it there).

However, a full swap of `ParSet`'s *base class* to `pydantic.BaseModel`
turns out not to be "thin": a codebase survey found dict-style
`par['a']['b']` get/set used across **~109 files** (thousands of call
sites, e.g. `pypeit/pypeit_steps.py`, `alignframe.py`, every
`pypeit/spectrographs/*.py`), `.parameters[key]['descr'/'options']`
metadata read directly in `pypeit/scripts/collate_1d.py`/`tellfit.py` and
by `to_rst_table()`, **15 of 41 subclasses** using the
`FrameGroupPar.parameters | {...}` dict-merge idiom (which would become
ordinary class inheritance under pydantic — an improvement, but a
mechanical edit to all 15), and a delicate ConfigObj `.pypeit`-text +
FITS-header I/O stack (`config_lines`/`to_config`, `recursive_dict_evaluate`
string coercion, `to_header`/`from_header`) that must be preserved
byte-for-byte regardless of what sits underneath it, and would need
porting either way. Altogether this is comparable in size to the rewrite
just completed on this branch, not a bolt-on.

**Decision (confirmed with user): pursue the small-footprint path only.**
Keep `ParSet` exactly as designed on this branch; apply the pydantic-derived
ideas (an explicit `default_factory` kwarg, `__init_subclass__` instead of a
metaclass for `FuncPar`) as isolated fixes without introducing
`pydantic.BaseModel` into `pypeitpar.py`. A full pydantic-backed rewrite (or
the pydantic-dataclasses hybrid considered as a middle ground) is
explicitly deferred, not rejected outright — worth revisiting only as its
own deliberate, separately-scoped decision, not as a follow-on to this
pass.

## Not in scope for this pass

- No wholesale replacement of `ParSet` with `attrs`/`pydantic`/`traitlets`
  (see feasibility exploration above).
- Porting `FuncPar` itself off the `funcpar_update`/`sensfunc_merge`
  branches into `parset_refactor-rebase` is separate follow-up work (noted
  in the standing plan for rebasing branches stacked on `parset_refactor`);
  this pass only proposes how `FuncPar` itself should be built once it is
  ported.

## First step on approval

Save this document verbatim to
`PypeIt-development-suite/claude_prompts/parset_refactor/alternate_approach.md`
(requested destination for this assessment) before making any code changes.

## Verification (once implemented)

- `pypeit/tests/test_pypeitpar.py` plus a small new direct test:
  instantiate `ReduxPar()` standalone (no `PypeItPar` wrapper) and assert
  `isinstance(par['redux_path'], (str, Path))` immediately after
  construction, with no `fill_callable()`/`validate()` call needed.
- Exercise `to_header()`/`to_config()` on a standalone `ReduxPar()` to
  confirm no callable/method repr leaks into the output.
- If `FuncPar` is ported in this pass: instantiate `DifferentialEvolutionPar()`
  and confirm `.parameters` matches `inspect.signature(scipy.optimize.differential_evolution)`
  minus `omitted_keys`, that `dtype` reflects annotated parameter types where
  available, and that no metaclass is involved (`type(DifferentialEvolutionPar)
  is type`).
