document.addEventListener('DOMContentLoaded', function () {
  const sections = document.querySelectorAll('.variable-section');

  sections.forEach(function (section) {
    const tabs = section.querySelectorAll('.tab');
    const panels = section.querySelectorAll('.panel');

    // Projection-specific links for this section only.
    const projectionLinks = section.querySelectorAll(
      '.code-link, .youtube-link'
    );

    tabs.forEach(function (tab) {
      tab.addEventListener('click', function () {
        const target = tab.getAttribute('data-projection');

        // Update this section's tabs only.
        tabs.forEach(function (t) {
          const isMatch = t === tab;

          t.classList.toggle('is-active', isMatch);
          t.setAttribute(
            'aria-selected',
            isMatch ? 'true' : 'false'
          );
        });

        // Show the selected projection panel.
        panels.forEach(function (panel) {
          const isMatch =
            panel.getAttribute('data-projection') === target;

          panel.classList.toggle('is-active', isMatch);

          // Pause videos in hidden panels.
          if (!isMatch) {
            const video = panel.querySelector('video');

            if (video && !video.paused) {
              video.pause();
            }
          }
        });

        // Update Python, R, and YouTube links.
        projectionLinks.forEach(function (link) {
          const url = link.getAttribute('data-' + target);

          if (url) {
            link.href = url;
          }
        });
      });
    });
  });

  // Playback speed controls.
  const speedButtons = document.querySelectorAll('.speed-btn');

  speedButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel
        ? panel.querySelector('video')
        : null;

      if (!video) return;

      video.playbackRate = parseFloat(
        btn.getAttribute('data-speed')
      );

      const siblingButtons =
        panel.querySelectorAll('.speed-btn');

      siblingButtons.forEach(function (b) {
        b.classList.toggle('is-active', b === btn);
      });
    });
  });

  // Loop controls.
  const loopButtons =
    document.querySelectorAll('[data-loop-toggle]');

  loopButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel
        ? panel.querySelector('video')
        : null;

      if (!video) return;

      video.loop = !video.loop;

      btn.classList.toggle(
        'is-active',
        video.loop
      );

      btn.setAttribute(
        'aria-pressed',
        video.loop ? 'true' : 'false'
      );
    });
  });

  // Read more / Show less controls.
  const aboutWraps =
    document.querySelectorAll('.meta-about-wrap');

  aboutWraps.forEach(function (wrap) {
    const text = wrap.querySelector('.meta-about');
    const toggle =
      wrap.querySelector('.meta-about-toggle');

    if (!text || !toggle) return;

    const isOverflowing =
      text.scrollHeight > text.clientHeight + 1;

    if (!isOverflowing) {
      toggle.style.display = 'none';
      return;
    }

    toggle.addEventListener('click', function () {
      const expanded =
        text.classList.toggle('is-expanded');

      toggle.textContent = expanded
        ? 'Show less'
        : 'Read more';
    });
  });
});
// Video download links: point each section's link at the active panel's video
document.querySelectorAll('.variable-section').forEach((section) => {
  const link = section.querySelector('.video-download');
  if (!link) return;

  const update = () => {
    const source = section.querySelector('.panel.is-active video source');
    if (source) link.href = source.getAttribute('src');
  };

  update();                                  // set on page load
  link.addEventListener('click', update);    // refresh at click time, after any tab change
});

// Sticky sidebar: built from the page's category headings and sections,
// with scroll-spy highlighting. New animations appear automatically as long
// as they use .category-heading and .variable-section (with .section-heading h2).
(function () {
  const nav = document.getElementById('side-nav');
  const toggle = document.getElementById('side-nav-toggle');
  if (!nav) return;

  const categories = [];
  const spyItems = []; // flat, document order: { el, group, link }
  const seen = new Set();
  let current = null;

  function startCategory(heading) {
    if (seen.has(heading)) return;
    seen.add(heading);
    const h2 = heading.querySelector('h2');
    current = {
      id: heading.id,
      title: h2 ? h2.textContent.trim() : heading.id,
      el: heading,
      sections: []
    };
    categories.push(current);
    spyItems.push({ el: heading, cat: current, sec: null });
  }

  document.querySelectorAll('.category-heading, .variable-section').forEach(function (el) {
    if (el.classList.contains('category-heading')) {
      startCategory(el);
      return;
    }
    // A section may carry its own category heading as a child (e.g. #lst).
    const inner = el.querySelector(':scope > .category-heading');
    if (inner) startCategory(inner);
    if (!current) return;
    const h2 = el.querySelector('.section-heading h2');
    if (!h2 || !el.id) return;
    const sec = { id: el.id, title: h2.textContent.replace(/\s+/g, ' ').trim(), el: el };
    current.sections.push(sec);
    spyItems.push({ el: el, cat: current, sec: sec });
  });

  // Build the markup.
  const title = document.createElement('a');
  title.className = 'side-nav-title';
  title.href = '#';
  title.textContent = 'Animation Library';
  nav.appendChild(title);

  const list = document.createElement('ul');
  categories.forEach(function (cat) {
    const li = document.createElement('li');
    li.className = 'side-nav-group';

    const a = document.createElement('a');
    a.className = 'side-nav-cat';
    a.href = '#' + cat.id;
    const label = document.createElement('span');
    label.textContent = cat.title;
    const count = document.createElement('span');
    count.className = 'side-nav-count';
    count.textContent = cat.sections.length;
    a.append(label, count);
    li.appendChild(a);

    const sub = document.createElement('ul');
    sub.className = 'side-nav-items';
    cat.sections.forEach(function (sec) {
      const sli = document.createElement('li');
      const sa = document.createElement('a');
      sa.className = 'side-nav-link';
      sa.href = '#' + sec.id;
      sa.textContent = sec.title;
      sec.link = sa;
      sli.appendChild(sa);
      sub.appendChild(sli);
    });
    li.appendChild(sub);
    cat.group = li;
    list.appendChild(li);
  });
  nav.appendChild(list);

  // Mobile / narrow: toggle drawer.
  function setOpen(open) {
    nav.classList.toggle('is-open', open);
    if (toggle) toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
  }
  if (toggle) {
    toggle.addEventListener('click', function () {
      setOpen(!nav.classList.contains('is-open'));
    });
  }
  nav.addEventListener('click', function (e) {
    if (e.target.closest('a')) setOpen(false);
  });
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') setOpen(false);
  });

  // Scroll-spy.
  const OFFSET = 160;
  let activeCat = null;
  let activeSec = null;
  let ticking = false;

  function update() {
    ticking = false;
    let hit = null;
    for (let i = 0; i < spyItems.length; i++) {
      if (spyItems[i].el.getBoundingClientRect().top <= OFFSET) hit = spyItems[i];
      else break;
    }
    const cat = hit ? hit.cat : null;
    const sec = hit ? hit.sec : null;
    if (cat === activeCat && sec === activeSec) return;

    if (activeSec && activeSec.link) {
      activeSec.link.classList.remove('is-active');
      activeSec.link.removeAttribute('aria-current');
    }
    if (activeCat && cat !== activeCat) {
      activeCat.group.classList.remove('is-active', 'is-open');
    }
    activeCat = cat;
    activeSec = sec;

    if (cat) cat.group.classList.add('is-active', 'is-open');
    if (sec && sec.link) {
      sec.link.classList.add('is-active');
      sec.link.setAttribute('aria-current', 'true');
    }

    // Keep the active entry visible inside the sidebar without moving the page.
    const target = (sec && sec.link) || (cat && cat.group.firstElementChild);
    if (target) {
      const n = nav.getBoundingClientRect();
      const t = target.getBoundingClientRect();
      if (t.top < n.top + 40) nav.scrollTop -= (n.top + 40 - t.top);
      else if (t.bottom > n.bottom - 40) nav.scrollTop += (t.bottom - n.bottom + 40);
    }
  }

  window.addEventListener('scroll', function () {
    if (!ticking) { ticking = true; requestAnimationFrame(update); }
  }, { passive: true });
  window.addEventListener('resize', update);
  update();
})();
