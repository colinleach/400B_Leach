--
-- PostgreSQL database dump
--

-- Dumped from database version 10.12 (Ubuntu 10.12-0ubuntu0.18.04.1)
-- Dumped by pg_dump version 10.12

-- Started on 2020-03-15 16:35:38

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- TOC entry 200 (class 1259 OID 70266)
-- Name: angmom; Type: TABLE; Schema: public; Owner: python
--

CREATE TABLE public.angmom (
    gal character(3) NOT NULL,
    snap smallint NOT NULL,
    t real NOT NULL,
    x_hat real NOT NULL,
    y_hat real NOT NULL,
    z_hat real NOT NULL,
    l_mag real NOT NULL
);


ALTER TABLE public.angmom OWNER TO python;

--
-- TOC entry 2923 (class 0 OID 0)
-- Dependencies: 200
-- Name: TABLE angmom; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON TABLE public.angmom IS 'Angular momentum of disk stars';


--
-- TOC entry 2924 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.gal; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.gal IS '''MW'', ''M31'', ''M33''';


--
-- TOC entry 2925 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.t; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.t IS 'Elapsed time, Myr';


--
-- TOC entry 2926 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.x_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.x_hat IS 'unit vector component';


--
-- TOC entry 2927 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.y_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.y_hat IS 'unit vector component';


--
-- TOC entry 2928 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.z_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.z_hat IS 'unit vector component';


--
-- TOC entry 2929 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.l_mag; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.l_mag IS 'magnitude of angular momentum';


--
-- TOC entry 2796 (class 2606 OID 70270)
-- Name: angmom angmom_pkey; Type: CONSTRAINT; Schema: public; Owner: python
--

ALTER TABLE ONLY public.angmom
    ADD CONSTRAINT angmom_pkey PRIMARY KEY (gal, snap);


-- Completed on 2020-03-15 16:35:39

--
-- PostgreSQL database dump complete
--

