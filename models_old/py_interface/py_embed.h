#ifndef PY_EMBED_CLASS_H
#define PY_EMBED_CLASS_H


class PyEmbed {
    PyEmbed() { py::initialize_interpreter(); }
    virtual ~PyEmbed() { py::finalize_interpreter(); }
};


#endif  // PY_EMBED_CLASS_H