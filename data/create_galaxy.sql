--
-- PostgreSQL database dump
--

-- Dumped from database version 10.12 (Ubuntu 10.12-0ubuntu0.18.04.1)
-- Dumped by pg_dump version 10.12

-- Started on 2020-03-15 16:37:53

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

--
-- TOC entry 1 (class 3079 OID 13039)
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- TOC entry 2936 (class 0 OID 0)
-- Dependencies: 1
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


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
-- TOC entry 2937 (class 0 OID 0)
-- Dependencies: 200
-- Name: TABLE angmom; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON TABLE public.angmom IS 'Angular momentum of disk stars';


--
-- TOC entry 2938 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.gal; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.gal IS '''MW'', ''M31'', ''M33''';


--
-- TOC entry 2939 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.t; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.t IS 'Elapsed time, Myr';


--
-- TOC entry 2940 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.x_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.x_hat IS 'unit vector component';


--
-- TOC entry 2941 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.y_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.y_hat IS 'unit vector component';


--
-- TOC entry 2942 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.z_hat; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.z_hat IS 'unit vector component';


--
-- TOC entry 2943 (class 0 OID 0)
-- Dependencies: 200
-- Name: COLUMN angmom.l_mag; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.angmom.l_mag IS 'magnitude of angular momentum';


--
-- TOC entry 198 (class 1259 OID 70208)
-- Name: centerofmass; Type: TABLE; Schema: public; Owner: python
--

CREATE TABLE public.centerofmass (
    snap smallint DEFAULT 0 NOT NULL,
    t real NOT NULL,
    x real NOT NULL,
    y real NOT NULL,
    z real NOT NULL,
    vx real NOT NULL,
    vy real NOT NULL,
    vz real NOT NULL,
    gal character(3) NOT NULL
);


ALTER TABLE public.centerofmass OWNER TO python;

--
-- TOC entry 2944 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.snap; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.snap IS 'primary key';


--
-- TOC entry 2945 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.t; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.t IS 'elapsed time, Myr';


--
-- TOC entry 2946 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.x; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.x IS 'position, kpc';


--
-- TOC entry 2947 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.y; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.y IS 'position, kpc';


--
-- TOC entry 2948 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.z; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.z IS 'position, kpc';


--
-- TOC entry 2949 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.vx; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.vx IS 'velocity, km/s';


--
-- TOC entry 2950 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.vy; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.vy IS 'velocity, km/s';


--
-- TOC entry 2951 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN centerofmass.vz; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.centerofmass.vz IS 'velocity, km/s';


--
-- TOC entry 199 (class 1259 OID 70211)
-- Name: centerofmass_inx_seq; Type: SEQUENCE; Schema: public; Owner: python
--

CREATE SEQUENCE public.centerofmass_inx_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.centerofmass_inx_seq OWNER TO python;

--
-- TOC entry 2952 (class 0 OID 0)
-- Dependencies: 199
-- Name: centerofmass_inx_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: python
--

ALTER SEQUENCE public.centerofmass_inx_seq OWNED BY public.centerofmass.snap;


--
-- TOC entry 196 (class 1259 OID 70111)
-- Name: simdata; Type: TABLE; Schema: public; Owner: python
--

CREATE TABLE public.simdata (
    inx integer NOT NULL,
    galname character(3) NOT NULL,
    snap smallint NOT NULL,
    "time" real NOT NULL,
    x real NOT NULL,
    y real NOT NULL,
    z real NOT NULL,
    vx real NOT NULL,
    vy real NOT NULL,
    vz real NOT NULL,
    m real NOT NULL,
    type smallint NOT NULL,
    pnum integer NOT NULL
);


ALTER TABLE public.simdata OWNER TO python;

--
-- TOC entry 2953 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.galname; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.galname IS 'MW, M31, M33';


--
-- TOC entry 2954 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.snap; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.snap IS '0 to 801';


--
-- TOC entry 2955 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata."time"; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata."time" IS 'units of Myr';


--
-- TOC entry 2956 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.x; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.x IS 'position, kpc';


--
-- TOC entry 2957 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.y; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.y IS 'position, kpc';


--
-- TOC entry 2958 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.z; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.z IS 'position, kpc';


--
-- TOC entry 2959 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.vx; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.vx IS 'velocity, km/s';


--
-- TOC entry 2960 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.vy; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.vy IS 'velocity, km/s';


--
-- TOC entry 2961 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.vz; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.vz IS 'velocity, km/s';


--
-- TOC entry 2962 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.m; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.m IS 'units of 10^10 M_sun';


--
-- TOC entry 2963 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.type; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.type IS '1: halo; 2: disk; 3: bulge';


--
-- TOC entry 2964 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN simdata.pnum; Type: COMMENT; Schema: public; Owner: python
--

COMMENT ON COLUMN public.simdata.pnum IS 'particle id (for creating a unique key)';


--
-- TOC entry 197 (class 1259 OID 70114)
-- Name: simdata_inx_seq; Type: SEQUENCE; Schema: public; Owner: python
--

CREATE SEQUENCE public.simdata_inx_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.simdata_inx_seq OWNER TO python;

--
-- TOC entry 2965 (class 0 OID 0)
-- Dependencies: 197
-- Name: simdata_inx_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: python
--

ALTER SEQUENCE public.simdata_inx_seq OWNED BY public.simdata.inx;


--
-- TOC entry 201 (class 1259 OID 70278)
-- Name: totalcom; Type: TABLE; Schema: public; Owner: python
--

CREATE TABLE public.totalcom (
    t real NOT NULL,
    x real NOT NULL,
    y real NOT NULL,
    z real NOT NULL,
    vx real NOT NULL,
    vy real NOT NULL,
    vz real NOT NULL,
    snap smallint
);


ALTER TABLE public.totalcom OWNER TO python;

--
-- TOC entry 2799 (class 2604 OID 70116)
-- Name: simdata inx; Type: DEFAULT; Schema: public; Owner: python
--

ALTER TABLE ONLY public.simdata ALTER COLUMN inx SET DEFAULT nextval('public.simdata_inx_seq'::regclass);


--
-- TOC entry 2807 (class 2606 OID 70270)
-- Name: angmom angmom_pkey; Type: CONSTRAINT; Schema: public; Owner: python
--

ALTER TABLE ONLY public.angmom
    ADD CONSTRAINT angmom_pkey PRIMARY KEY (gal, snap);


--
-- TOC entry 2805 (class 2606 OID 70277)
-- Name: centerofmass primary; Type: CONSTRAINT; Schema: public; Owner: python
--

ALTER TABLE ONLY public.centerofmass
    ADD CONSTRAINT "primary" PRIMARY KEY (gal, snap);


--
-- TOC entry 2803 (class 2606 OID 70119)
-- Name: simdata simdata_pkey; Type: CONSTRAINT; Schema: public; Owner: python
--

ALTER TABLE ONLY public.simdata
    ADD CONSTRAINT simdata_pkey PRIMARY KEY (galname, snap, pnum);


--
-- TOC entry 2801 (class 1259 OID 70120)
-- Name: simdata_inx; Type: INDEX; Schema: public; Owner: python
--

CREATE UNIQUE INDEX simdata_inx ON public.simdata USING btree (inx);


-- Completed on 2020-03-15 16:37:53

--
-- PostgreSQL database dump complete
--

