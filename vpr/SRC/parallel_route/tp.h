#undef TRACEPOINT_PROVIDER
#define TRACEPOINT_PROVIDER hello_world

#undef TRACEPOINT_INCLUDE
#define TRACEPOINT_INCLUDE "./tp.h"

#if !defined(_HELLO_TP_H) || defined(TRACEPOINT_HEADER_MULTI_READ)
#define _HELLO_TP_H

#include <lttng/tracepoint.h>

TRACEPOINT_EVENT(
    hello_world,
    my_first_tracepoint,
    TP_ARGS(
        int, tid,
		int, context,
		int, type,
		int, lock_count
    ),
    TP_FIELDS(
        ctf_integer(int, tid_field, tid)
        ctf_integer(int, context_field, context)
        ctf_integer(int, type_field, type)
        ctf_integer(int, lock_count_field, lock_count)
    )
)

#endif /* _HELLO_TP_H */

#include <lttng/tracepoint-event.h>
